#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Moving sky plot to bottom of object spectrum plot, added badlines
# stuff,restricted z range
# for cross correlation to slider z +/-0.3
#
# Adding xcsao button - need shady iraf import shenanigans, otherwise iraf
# causes feedback between tk
# widgets and variables to not work!
#
# Same as visualTemplateRedshift.py, but accepts wildcards and automatically
# saves out plots/results
# (quit button -> done button)
#
# Has a lovely tk interface
#
# Overplots a spectrum and spectral template using astSED functions
# This works with SDSS log wavelength templates
#
# Input spectra can either be DEEP2 pipeline spec1d files, or similar
# formatted FITS tables with
# spectrum in extension '1D_SPECTRUM'. Masking out values in input spectrum
# with values > 10*median
# (gets rid of dodgy ends in efosc2 spectra)

import os
import sys
import datetime
import Tkinter
import pyfits
import numpy
import pylab
import matplotlib.patches as patches
from astLib import astSED
from astLib import astStats
from astLib import astWCS
from scipy import optimize
from scipy import interpolate

pylab.matplotlib.interactive(True)

#-----------------------------------------------------------------------------
# SDSS templates to use - note later on we explicitly load in the Tremonti et
# al. starburst template
# which is in a different format and stored under
# $HOME/Astro_Software/TremontiStarburstTemplate
tempDir = "./VisualTemplateRedshiftTemplates"
#tremontiFileName=os.environ['HOME']+os.path.sep+"Astro_Software"+os.path.sep+"TremontiStarburstTemplate/fos_ghrs_composite.txt"

templateFileNames = [
    # Galaxies
    "spDR2-023.fit",
    "spDR2-024.fit",
    "spDR2-025.fit",
    "spDR2-026.fit",
    "spDR2-027.fit",
    "spDR2-028.fit",
    # QSOs
    "spDR2-029.fit",
    "spDR2-030.fit",
    "spDR2-031.fit",
    "spDR2-032.fit"
]
# Added this because Tremonti starburst template added
templateLabels = [
    # Galaxies
    "SDSS-023",
    "SDSS-024",
    "SDSS-025",
    "SDSS-026",
    "SDSS-027",
    "SDSS-028",
    # QSOs
    "SDSS-029",
    "SDSS-030",
    "SDSS-031",
    "SDSS-032"
    #,
    # Starburst
    #"T03 Starburst"
]

# List of spectral line labels and wavelengths
# Unlike the old style scripts, these don't need unique labels/names
spectralFeaturesCatalogue = [
    ["Halpha", 6562.8],
    ["Hbeta", 4861.33],
    ["Hgamma", 4340.50],
    ["Hdelta", 4101.70],
    ["Hepsilon", 3970.10],
    ["Htheta", 3798.6],
    ["Hzeta", 3888.9],
    ["HeI", 4471.5, 5875.6],
    ["HeII", 4338.6, 4685.70],
    ["Lyalpha", 1216.0],
    ["CII", 1335.30],
    ["CIV", 1549.0],
    ["[NII]", 6548.1, 6583.4],
    ["NV", 1240.14],
    ["OI", 1304.3],
    ["[OI]", 6300.3],
    ["[OII]", 2471.03, 3727.3],
    ["OIII", 2672.04],
    ["[OIII]", 4958.91, 5006.84],
    ["OIV]", 1402.06],
    ["[NeIII]", 3967.5],
    ["[NeIV]", 1602.0, 2423.83],
    ["[NeV]", 1575.0, 3426.0],
    ["MgI", 2852.0, 3830.4, 5175.4],
    ["MgII", 2800.0, 2803.0, 2798.0],
    ["AlII]", 2669.95],
    ["SiII", 1262.59, 1306.82],
    ["SiIV", 1396.76],
    ["[SII]", 6717.0, 6731.3],
    ["CaI", 4226.7],
    ["H", 3968.47],
    ["K", 3933.68],
    ["G", 4307.74],
    ["E", 5269.00],
    ["FeI", 4045.8, 4271.7, 4329.7, 4045.8],
    ["FeII", 2344.0, 2374.0, 2382.0, 2586.0, 2600.0]
]


class App:
    def __init__(self, master, objectSpecFileNames, outDir):
        """outDir is a dir in which to write output

        """
        # Do some initial setup of things
        self.objectSpecFileNames = objectSpecFileNames
        self.currentSpecFileIndex = 0

        # Create the output text file for results
        self.outDir = outDir
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        #self.outFile = file(outDir + os.path.sep +
        #                    datetime.datetime.today().isoformat() +
        #                    ".results", "w")
        self.outFile = file(outDir + os.path.sep +\
            self.objectSpecFileNames[self.currentSpecFileIndex].rstrip('.fits')
            + ".results", "w")
        self.outFile.write("#Fiber\tRedshift\tRedshiftError" +
                           "\tQuality\tComments\n")
        self.outFile.close()
        # Load Templates
        self.templates = self.loadTemplates(tempDir, templateFileNames)

        # Load Science
        print "Checking spectrum %d/%d ..." % (self.currentSpecFileIndex + 1,
                                               len(self.objectSpecFileNames))
        # Single Frame
        if 0:
            obj = self.loadObjectSpectrum(objectSpecFileNames[
                self.currentSpecFileIndex])
            self.objectSpecFileName = objectSpecFileNames[
                self.currentSpecFileIndex]  # for plot titles etc.
            self.objSED = obj['object']
            self.unsmoothedObjFlux = self.objSED.flux[:]
            self.skySED = obj['sky']
        # IFU data
        else:
            self.obj = self.loadIFUSpectra(objectSpecFileNames[
                self.currentSpecFileIndex])
            self.objectSpecFileName = objectSpecFileNames[
                self.currentSpecFileIndex]  # for plot titles etc.

        # List of feature names in spectralFeaturesCatalogue to plot
        self.plotFeatures = []

        ###########################
        ### BEGIN LAYOUT OF APP ###
        ###########################

        # Layout - these get drawn top to bottom
        self.buttonFrame = Tkinter.Frame(master, padx=5, pady=5)
        self.buttonFrame.grid()

        self.qualityFrame = Tkinter.Frame(master, padx=5, pady=5)
        self.qualityFrame.grid()

        self.fiberFrame = Tkinter.Frame(master, padx=5, pady=5)
        self.fiberFrame.grid()

        self.smoothFrame = Tkinter.Frame(master, padx=5, pady=5)
        self.smoothFrame.grid()

        templatesFrame = Tkinter.Frame(master, padx=5, pady=5)
        templatesFrame.grid()

        scaleFrame = Tkinter.Frame(master, padx=5, pady=5)
        scaleFrame.grid()

        featuresFrame = Tkinter.Frame(master, padx=5, pady=5)
        featuresFrame.grid()

        # Buttons frame
        self.redrawButton = Tkinter.Button(self.buttonFrame,
                                           text="Redraw plot",
                                           command=self.redrawPlot)
        self.redrawButton.grid(row=0, column=0)

        self.savePNGButton = Tkinter.Button(self.buttonFrame,
                                            text="Save .png",
                                            command=self.savePNG)
        self.savePNGButton.grid(row=0, column=1)

        self.outPathLabel = Tkinter.Label(self.buttonFrame,
                                          text="Output .png file : ")
        self.outPathLabel.grid(row=0, column=2)
        self.outPathEntryVar = Tkinter.StringVar()
        self.outPathEntryVar.set(self.outDir + os.path.sep +
                                 self.objectSpecFileName.replace(
                                     ".fits", ".png"))
        self.outPathEntry = Tkinter.Entry(self.buttonFrame,
                                          textvariable=self.outPathEntryVar,
                                          width=80)
        self.outPathEntry.grid(row=0, column=3)

        self.nextButton = Tkinter.Button(self.buttonFrame,
                                         text="Log",
                                         bg="blue", command=self.log)
        self.nextButton.grid(row=0, column=6)

        self.quitButton = Tkinter.Button(self.buttonFrame,
                                         text="QUIT", bg="red",
                                         command=self.buttonFrame.quit)
        self.quitButton.grid(row=0, column=7)

        # Quality frame, contains radio buttons and comments field
        self.qualityRadioVar = Tkinter.IntVar()
        self.qualityRadioList = []
        self.qualityLabel = Tkinter.Label(self.qualityFrame,
                                          text="Quality flag : ",
                                          anchor=Tkinter.E)
        self.qualityLabel.grid(row=2, column=0)
        for i in range(4):
            self.qualityRadioList.append(
                Tkinter.Radiobutton(self.qualityFrame,
                                    text=str(i),
                                    variable=self.qualityRadioVar,
                                    value=i,
                                    command=self.redrawPlot))
            self.qualityRadioList[-1].grid(row=2, column=i + 1)
        self.qualityRadioList[0].select()

        # Comment Box
        self.commentsLabel = Tkinter.Label(self.qualityFrame,
                                           text="Comments: ")
        self.commentsLabel.grid(row=2, column=i + 2)
        self.commentsEntryVar = Tkinter.StringVar()
        self.commentsEntryVar.set("")
        self.commentsEntry = Tkinter.Entry(self.qualityFrame,
                                           textvariable=self.commentsEntryVar,
                                           width=80)
        self.commentsEntry.grid(row=2, column=i + 3)

        # Slider to choose fiber
        self.fibernumberVar = Tkinter.DoubleVar()
        self.fibernumberLabel = Tkinter.Label(self.fiberFrame,
                                              text='Fiber Number : ',
                                              anchor=Tkinter.E)
        self.fibernumberLabel.grid(row=1, column=0)
        self.fibernumberScale = Tkinter.Scale(self.fiberFrame,
                                              orient=Tkinter.HORIZONTAL,
                                              length=300, from_=1, to=246,
                                              tickinterval=60,
                                              command=None,
                                              variable=self.fibernumberVar,
                                              resolution=1)
        self.fibernumberScale.set(1)
        self.fibernumberScale.grid(row=1, column=1, columnspan=3)

        # Buttons to finely tune the fiber number
        self.fibernumberMinusButton = Tkinter.Button(self.fiberFrame,
                                                     text="-",
                                                     command=self.decreaseFiber)
        self.fibernumberMinusButton.grid(row=1, column=4)
        self.fibernumberPlusButton = Tkinter.Button(self.fiberFrame,
                                                    text="+",
                                                    command=self.increaseFiber)
        self.fibernumberPlusButton.grid(row=1, column=5)

        # Checkbox to set fiber as sky
        self.ignoreEmission = Tkinter.IntVar()
        self.ignoreEmissionLabel = Tkinter.Label(self.fiberFrame,
                            text="Ignore Emission", width=15, anchor=Tkinter.E)
        self.ignoreEmissionLabel.grid(row=1, column=10)
        self.ignoreEmissionCheckButton = Tkinter.Checkbutton(
            self.fiberFrame,
            variable=self.ignoreEmission,
            command=None)
        self.ignoreEmissionCheckButton.grid(row=1, column=11)

        # Slider to set smoothing of object spectrum
        self.smoothScaleVar = Tkinter.DoubleVar()
        self.smoothLabel = Tkinter.Label(self.smoothFrame,
                                         text="Spectrum smoothing : ",
                                         anchor=Tkinter.E)
        self.smoothLabel.grid(row=1, column=0)
        self.smoothScale = Tkinter.Scale(self.smoothFrame,
                                         orient=Tkinter.HORIZONTAL,
                                         length=300, from_=0, to=100,
                                         tickinterval=25,
                                         command=self.smoothSpectrum,
                                         variable=self.smoothScaleVar,
                                         resolution=1)

        # Initial Smoothing
        self.smoothScale.set(5)
        x = self.smoothScale.get()
        self.objSED = self.obj[246-self.fibernumberScale.get()]['object']
        self.unsmoothedObjFlux = self.objSED.flux[:]
        self.objSED.flux = self.unsmoothedObjFlux[:]
        self.objSED.smooth(x)
        self.smoothScale.grid(row=1, column=1, columnspan=3)

        # Buttons to finely tune the smoothing
        self.smoothPlusButton = Tkinter.Button(self.smoothFrame, text="+",
                                               command=self.increaseSmoothing)
        self.smoothPlusButton.grid(row=1, column=5)
        self.smoothMinusButton = Tkinter.Button(self.smoothFrame, text="-",
                                                command=self.decreaseSmoothing)
        self.smoothMinusButton.grid(row=1, column=4)

        # Min, max wavelength range
        self.minWavelengthLabel = Tkinter.Label(self.smoothFrame,
                                                text="Min WL:",
                                                anchor=Tkinter.E)
        self.minWavelengthLabel.grid(row=1, column=6)
        self.minWavelengthEntryVar = Tkinter.StringVar()
        self.minWavelengthEntryVar.set("4500")
        self.minWavelengthEntry = Tkinter.Entry(
            self.smoothFrame,
            textvariable=self.minWavelengthEntryVar, width=6)
        self.minWavelengthEntry.grid(row=1, column=7)

        self.maxWavelengthLabel = Tkinter.Label(self.smoothFrame,
                                                text="Max WL:",
                                                anchor=Tkinter.E)
        self.maxWavelengthLabel.grid(row=1, column=8)
        self.maxWavelengthEntryVar = Tkinter.StringVar()
        self.maxWavelengthEntryVar.set("6500")
        self.maxWavelengthEntry = Tkinter.Entry(
            self.smoothFrame,
            textvariable=self.maxWavelengthEntryVar, width=6)
        self.maxWavelengthEntry.grid(row=1, column=9)

        # Alt normalisation method
        self.altNormCheckVar = Tkinter.IntVar()
        self.altNormLabel = Tkinter.Label(self.smoothFrame, text="Alt norm",
                                          width=10, anchor=Tkinter.E)
        self.altNormLabel.grid(row=1, column=10)
        self.altNormCheckButton = Tkinter.Checkbutton(
            self.smoothFrame,
            variable=self.altNormCheckVar, command=self.plotSkyChanged)
        self.altNormCheckButton.grid(row=1, column=11)

        # Turn sky plotting on/off
        self.plotSkyCheckVar = Tkinter.IntVar()
        self.plotSkyLabel = Tkinter.Label(self.smoothFrame, text="Plot sky",
                                          width=10, anchor=Tkinter.E)
        self.plotSkyLabel.grid(row=1, column=12)
        self.plotSkyCheckButton = Tkinter.Checkbutton(
            self.smoothFrame,
            variable=self.plotSkyCheckVar,
            command=self.plotSkyChanged)
        self.plotSkyCheckButton.grid(row=1, column=13)

        # Templates frame
        # Radio buttons used to select template
        self.templateRadioVar = Tkinter.IntVar()
        self.templateRadioList = []
        self.templateLabel = Tkinter.Label(templatesFrame, text="Template:",
                                           anchor=Tkinter.E)
        self.templateLabel.grid(row=2, column=0)

        for i in range(len(self.templates)):
            try:
                tempName = templateLabels[i]
            except IndexError:
                sys.exit()
            self.templateRadioList.append(Tkinter.Radiobutton(
                templatesFrame,
                text=tempName,
                variable=self.templateRadioVar,
                value=i,
                command=self.redrawPlot))
            self.templateRadioList[-1].grid(row=2, column=i + 1)

        self.templateRadioList[0].select()

        # Slider used to set trial redshift of template
        self.redshiftScaleVar = Tkinter.DoubleVar()
        self.redshiftLabel = Tkinter.Label(scaleFrame,
                                           text="Template redshift : ",
                                           anchor=Tkinter.E)
        self.redshiftLabel.grid(row=3, column=0)
        self.redshiftScale = Tkinter.Scale(scaleFrame,
                                           orient=Tkinter.HORIZONTAL,
                                           length=600, from_=0.0, to=2.01,
                                           tickinterval=1,
                                           command=self.getRedshiftScaleValue,
                                           variable=self.redshiftScaleVar,
                                           resolution=0.001)
        self.redshiftScale.set(0.5)
        self.redshiftScale.grid(row=3, column=1, columnspan=3)

        # Buttons to finely tune the trial redshift
        self.redshiftPlusButton = Tkinter.Button(scaleFrame, text="+",
                                                 command=self.increaseRedshift)
        self.redshiftPlusButton.grid(row=3, column=5)
        self.redshiftMinusButton = Tkinter.Button(scaleFrame, text="-",
                                                  command=self.decreaseRedshift)
        self.redshiftMinusButton.grid(row=3, column=4)

        # Redshift uncertainty entry box
        self.redshiftErrorLabel = Tkinter.Label(scaleFrame, text="+/-")
        self.redshiftErrorLabel.grid(row=3, column=6)
        self.redshiftErrorEntryVar = Tkinter.StringVar()
        self.redshiftErrorEntryVar.set(0.001)
        self.redshiftErrorEntry = Tkinter.Entry(
            scaleFrame,
            textvariable=self.redshiftErrorEntryVar,
            width=10)
        self.redshiftErrorEntry.grid(row=3, column=7)

        # XCSAO button
        self.runXCSAOButton = Tkinter.Button(scaleFrame, text="XC galaxies",
                                             command=self.runXCSAOGalaxies,
                                             fg="green")
        self.runXCSAOButton.grid(row=3, column=8)
        self.runXCSAOButton = Tkinter.Button(scaleFrame, text="XC QSOs",
                                             command=self.runXCSAOQSOs,
                                             fg="green")
        self.runXCSAOButton.grid(row=3, column=9)

        # Features to optionally plot
        # Have to have separate labels for the check boxes because other
        # layout goes stupid
        self.featuresCheckList = []
        self.featuresCheckLabelsList = []
        self.featuresCheckVars = []
        maxPerRow = 20
        row = 0
        column = 0
        for i in range(len(spectralFeaturesCatalogue)):
            if column == maxPerRow:
                row += 1
                column = 0
            self.featuresCheckVars.append(Tkinter.IntVar())
            self.featuresCheckLabelsList.append(
                Tkinter.Label(featuresFrame,
                              text=spectralFeaturesCatalogue[i][0], width=10,
                              anchor=Tkinter.E))
            self.featuresCheckLabelsList[-1].grid(row=row + 4, column=column)
            self.featuresCheckList.append(Tkinter.Checkbutton(
                featuresFrame,
                variable=self.featuresCheckVars[-1],
                command=self.setPlotFeatures))
            self.featuresCheckList[-1].grid(row=row + 4, column=column + 1)
            column += 2

        # Start up the figure for drawing
        pylab.figure(figsize=(12, 8))

        # Do initial plot
        self.updatePlot(self.obj[246-self.fibernumberScale.get()]['object'],
                        self.templates[self.templateRadioVar.get()],
                        self.obj[246-self.fibernumberScale.get()]['sky'],
                        self.redshiftScaleVar.get(),
                        tempLabel=os.path.split(templateLabels[
                            self.templateRadioVar.get()])[-1],
                        redrawSky=True,
                        redrawFeatures=True,
                        plotFeatures=self.plotFeatures)

#        self.updatePlot(self.objSED,
#                        self.templates[self.templateRadioVar.get()],
#                        self.skySED,
#                        self.redshiftScaleVar.get(),
#                        tempLabel=os.path.split(templateLabels[
#                            self.templateRadioVar.get()])[-1],
#                        redrawSky=True,
#                        redrawFeatures=True,
#                        plotFeatures=self.plotFeatures)

    #########################
    ### END OF APP LAYOUT ###
    #########################
    def smoothSpectrum(self, event):
        """ Smooths the object spectrum when the smooth slider is updated

        """
        self.objSED = self.obj[246-self.fibernumberScale.get()]['object']
        #self.unsmoothedObjFlux = self.objSED.flux[:]
        #self.objSED.flux = self.unsmoothedObjFlux[:]
        self.objSED.smooth(self.smoothScale.get())

    def increaseSmoothing(self):
        """ Increases the smoothing by 1

        """
        self.smoothScale.set(self.smoothScaleVar.get() + 1)

    def decreaseSmoothing(self):
        """ Decreases the smoothing by 1

        """
        self.smoothScale.set(self.smoothScaleVar.get() - 1)

    def changeFiber(self):
        self.updatePlot(self.obj[246-self.fibernumberScale.get()]['object'],
                        self.templates[self.templateRadioVar.get()],
                        self.obj[246-self.fibernumberScale.get()]['sky'],
                        self.redshiftScaleVar.get(),
                        tempLabel=os.path.split(templateLabels[
                            self.templateRadioVar.get()])[-1],
                        redrawSky=True,
                        redrawFeatures=True,
                        plotFeatures=self.plotFeatures)

    def increaseFiber(self):
        """  Increases the fiber number by 1

        """
        self.fibernumberScale.set(self.fibernumberVar.get() + 1)
        self.outPathEntryVar.set(self.outDir + os.path.sep +
                                 self.objectSpecFileName.rstrip(".fits") +
                                 '_' + str(self.fibernumberScale.get()) +
                                 '.png')
        self.changeFiber()

    def decreaseFiber(self):
        """ Decreases the fiber number by 1

        """
        self.fibernumberScale.set(self.fibernumberVar.get() - 1)
        self.outPathEntryVar.set(self.outDir + os.path.sep +
                                 self.objectSpecFileName.rstrip(".fits") +
                                 '_' + str(self.fibernumberScale.get()) +
                                 '.png')
        self.changeFiber()

    def resetFeatures(self):
        """ Resets (turns off) the plotting of the spectral features

        """
        for i in range(len(self.plotFeatures)):
            self.featuresCheckVars[i].set(0)
        self.setPlotFeatures()

    def getRedshiftScaleValue(self, event):
        """ Gets the current value of the template redshift from the slider

        """
        self.redshiftScale.get()

    def increaseRedshift(self):
        """ Increases the template redshift by delta-z

        """
        self.redshiftScale.set(self.redshiftScaleVar.get() + 0.001)

    def decreaseRedshift(self):
        """ Decreases the template redshift by delta-z

        """
        self.redshiftScale.set(self.redshiftScaleVar.get() - 0.001)

    def setPlotFeatures(self):
        """ Sets self.plotFeatues and triggers redrawing of plot, according to
        which spectral features are selected. Redraws plot afterwards.

        """
        self.plotFeatures = []
        for i in range(len(self.featuresCheckVars)):
            val = self.featuresCheckVars[i].get()
            self.plotFeatures.append(val)
        self.redrawPlot()

    def savePNG(self):
        """ Saves the current figure to a png

        """
        pylab.savefig(self.outPathEntry.get())

    def runXCSAOGalaxies(self):
        """ Function which chooses what templates to include when matching
        galaxy type objects. Calls runXCSAO.

        """
        templatesToInclude = []
        for temp in templateFileNames:
            tempNum = int(os.path.split(temp)[-1].split("-")[-1].split(".")[0])
            if 23 <= tempNum < 29:
                templatesToInclude.append(temp)
        self.runXCSAO(templatesToInclude=templatesToInclude)

    def runXCSAOQSOs(self):
        """ Function which chooses what templates to include when matching
        QSO type objects. Calls runXCSAO.

        """
        templatesToInclude = []
        for temp in templateFileNames:
            tempNum = int(os.path.split(temp)[-1].split("-")[-1].split(".")[0])
            if 29 <= tempNum < 33:
                templatesToInclude.append(temp)
        self.runXCSAO(templatesToInclude=templatesToInclude)

    # This does all the work of fitting the spectra
    def runXCSAO(self, templatesToInclude=[]):
        try:
            from pyraf import iraf
        except ImportError:
            print 'There is something wrong with pyraf! Exiting...'
            sys.exit(1)
        try:
            from iraf import rvsao
        except ImportError:
            print 'RVSAO not installed? Exiting...'
            sys.exit(1)

        xMin = float(self.minWavelengthEntry.get())
        xMax = float(self.maxWavelengthEntry.get())

        # This has to have the file available in IRAF friendly format
        if not os.path.exists("spec1d_IRAF"):
            os.makedirs("spec1d_IRAF")
        irafFileName = "spec1d_IRAF" + os.path.sep + "iraf_boxcar_" + \
                      str(self.fibernumberScale.get())+ "_" +\
                      self.objectFileName
        print irafFileName
        result = None
        if os.path.exists(irafFileName) or \
                self.convertToIRAFFormat():
            print "--> cross correlating " + irafFileName + " ..."

            self.skySED = self.obj[246-self.fibernumberScale.get()]['sky']

            # Mask prominent sky emission lines
            if self.skySED is None:
                fixbad = "n"
            else:
                fixbad = "y"
            if fixbad == "y":
                normSkyFlux = self.skySED.flux / self.skySED.flux.max()
                threshold = 0.25
                badPix = numpy.where(normSkyFlux > threshold)[0]
                lines = []
                for i in range(len(badPix)):
                    startPixInLine = False
                    for line in lines:
                        if line[0] <= badPix[i] <= line[1]:
                            startPixInLine = True
                    if not startPixInLine:
                        pixCount = 1
                        if pixCount + i < len(badPix) - 1:
                            nextPix = badPix[i + pixCount]
                            prevPix = badPix[i]
                            maxReached = False
                            while nextPix < prevPix + 2:
                                if pixCount + i < len(badPix) - 1:
                                    prevPix = badPix[i + pixCount]
                                    pixCount += 1
                                    nextPix = badPix[i + pixCount]
                                else:
                                    maxReached = True
                                    break
                            if not maxReached:
                                lastPix = prevPix
                            else:
                                lastPix = max(badPix)
                        else:
                            lastPix = max(badPix)
                        lines.append([badPix[i], lastPix])
                    #outFile=file("/opt/iraf_extras/rvsao-2.8.1/lib/badlines.dat",
                #        "w")
                outFile = file("badlines.dat", "w")
                for line in lines:
                    print('Badlines updated')
                    outFile.write(str(self.skySED.wavelength[line[0]]-10) + \
                                  "\t" +\
                                  str(self.skySED.wavelength[line[1]]+10) +\
                                  "\n")
                outFile.close()

                # Cross correlate with SDSS galaxy templates
                #for zStep in range(0,9):
                #z=0.2*zStep

            if self.ignoreEmission.get() == 1:
                chop = 'y'
            else:
                chop = 'n'

            for temp in templatesToInclude:
                try:
                    rvsao.xcsao(spectra=irafFileName, tempdir=tempDir,
                                fixbad=fixbad,
                                badlines=os.getcwd()+os.path.sep+"badlines.dat",
                                vel_init="zguess",
                                czguess=self.redshiftScaleVar.get(),
                                templates=os.path.split(temp)[-1],
                                st_lambda=xMin,
                                end_lambda=xMax, zeropad="y", nsmooth=30,
                                #s_emchop="n", t_emchop="n", nzpass=8,
                                s_emchop=chop, t_emchop=chop, nzpass=8,
                                minvel=(self.redshiftScaleVar.get() - 0.3) *\
                                3e5,
                                maxvel=(self.redshiftScaleVar.get() + 0.3) *\
                                3e5,
                                renormalize="y", ncols=8192, low_bin=10,
                                top_low=20, top_nrun=250, nrun=500,
                                bell_window=0.05, dispmode=2, curmode="n",
                                ablines="ablines.dat", displot="no",
                                logfiles="xcsao.log", save_vel="n", pkfrac=0.5,
                                report_mode=1)
                except:
                    print "hmm? xcsao fell over"
                    sys.exit(1)
            result = self.parseXCSAOResult()
            os.remove("xcsao.log")

        if result is not None:
            self.redshiftScaleVar.set(result['z'])
            self.commentsEntryVar.set("XCSAO (R=%.3f): " % (result['RVal']))
            self.commentsEntry['textvariable'] = self.commentsEntryVar
            self.templateRadioVar.set(
                templateFileNames.index(result['template']))
            self.templateRadioList[
                templateFileNames.index(result['template'])].select()
            self.commentsEntry.update()
            self.redshiftErrorEntryVar.set("%.6f" % (result['zErr']))
            print result
        else:
            print "XCSAO failed."

        self.redrawPlot()

    def convertToIRAFFormat(self):
        """ Canned convert self.objectFileName spec1d IDL file into IRAF format,
        stored under 'spec1d_IRAF' dir. Needed for running XCSAO. Note if
        values < 1 here, we multiply by ridiculous factor to
        stop xcsao from crashing

        """
        try:
            from pyraf import iraf
        except ImportError:
            print 'There is something wrong with pyraf! Exiting...'
            sys.exit(1)
        try:
            from iraf import onedspec
        except ImportError:
            print 'Did not find onedspec in iraf! Exiting...'
            sys.exit(1)

        method = "boxcar"
        baseName = "spec1d_IRAF" + os.path.sep + "iraf_" + method + "_" \
                    + str(self.fibernumberScale.get()) + "_"
        print "--> Converting " + self.objectFileName + " to IRAF format ..."

        if method == "boxcar":
            tabExts = [1, 2]
        elif method == "optimal":
            tabExts = [3, 4]
        outFileName = baseName + self.objectFileName.replace(".fits", ".csv")
        writer = file(outFileName, "wb")
        skyWriter = file(outFileName.replace("iraf_", "sky_iraf_"), "wb")

        #idlfits = pyfits.open(self.objectFileName)

        fluxData = self.obj[246-self.fibernumberScale.get()]['object'].flux
        wavelengthData =\
        self.obj[246-self.fibernumberScale.get()]['object'].wavelength
        skyData = self.obj[246-self.fibernumberScale.get()]['sky'].flux

        print 246-self.fibernumberScale.get()
        print fluxData.mean()

        #if fluxData.mean() < 100:
        #    fluxData = (fluxData / fluxData.mean()) * 1e6
        plotData = []
        skyPlotData = []
        for i in range(len(fluxData)):
            writer.write(str(
                wavelengthData[i]) + "\t" + str(fluxData[i]) + "\n")
            skyWriter.write(str(
                wavelengthData[i])+"\t"+str(skyData[i])+"\n")
            plotData.append([wavelengthData[i], fluxData[i]])
            skyPlotData.append([wavelengthData[i], skyData[i]])
        writer.close()
        skyWriter.close()

        #for tabExt in tabExts:
        #    # Sometimes things just don't work as they should ...
        #    dataOkay = True
        #    try:
        #        if len(idlfits[tabExt].data.field('SPEC').shape) > 1:
        #            fluxData = idlfits[tabExt].data.field('SPEC')[0]
        #        else:
        #            fluxData = idlfits[tabExt].data.field('SPEC')
        #    except IndexError:
        #        dataOkay = False
        #    if dataOkay:
        #        if len(idlfits[tabExt].data.field('LAMBDA').shape) > 1:
        #            wavelengthData = idlfits[tabExt].data.field('LAMBDA')[0]
        #        else:
        #            wavelengthData = idlfits[tabExt].data.field('LAMBDA')
        #            #skyData=idlfits[tabExt].data.field('SKYSPEC')[0]
        #
        #        # xcsao crash preventing when running on fluxed spectra (it
        #        # should do this itself, of course)
        #        if fluxData.mean() < 100:
        #            fluxData = (fluxData / fluxData.mean()) * 1e6
        #        minWavelength = min(wavelengthData)
        #        maxWavelength = max(wavelengthData)
        #        plotData = []
        #        skyPlotData = []
        #        for i in range(len(fluxData)):
        #            writer.write(str(
        #                wavelengthData[i]) + "\t" + str(fluxData[i]) + "\n")
        #            #skyWriter.write(str(
        #            #    wavelengthData[i])+"\t"+str(skyData[i])+"\n")
        #            plotData.append([wavelengthData[i], fluxData[i]])
        #            #skyPlotData.append([wavelengthData[i],skyData[i]])
        #writer.close()
        #skyWriter.close()
        if not os.path.exists(outFileName.replace(".csv", ".fits")):
            #os.remove(outFile.rstrip("csv")+"fits")
            # This might fall over below if there's a dodgy file
            try:
                onedspec.rspectext(input=outFileName,
                                   output=outFileName.replace(".csv", ".fits"),
                                   flux="no", dtype="interp")
                onedspec.rspectext(input=outFileName.replace("iraf_",
                        "sky_iraf_"), output=outFileName.replace("iraf_",
                        "sky_iraf_").rstrip("csv")+"fits", flux="no",
                        dtype="interp")
                #del sys.modules['pyraf']
                #del sys.modules['pyraf.iraf']
                return True
            except:
                print "... Argh! there's a problem with this file that ", \
                    "causes IRAF to crash! ..."
                print "... skipping ..."

                return False

    def parseXCSAOResult(self):
        """ Parses xcsao results log file, returns highest R value result

        """
        inFile = file("xcsao.log", "rb")
        lines = inFile.readlines()
        #results = []
        objectList = []
        currentObject = ""
        bestRVal = 0.0
        bestResult = None
        for line in lines:
            # Check for spectrum name change
            objectChanged = False
            if line.find("Object:") != -1:
                newObject = line[line.find("Object:") + 8:].rstrip(" \n")
                if newObject != currentObject:
                    currentObject = newObject
                    objectChanged = True
            if objectChanged:
                objectList.append(currentObject)
                # Extract results -- if we asked to ignore any templates it's done
            # here
            if line.find("CZ:") != -1:
                bits = line.split()
                result = {}
                for i in range(len(bits)):
                    if bits[i] == "Temp:":
                        result['template'] = str(bits[i + 1])
                    if bits[i] == "R:":
                        result['RVal'] = float(bits[i + 1])
                    if bits[i] == "CZ:":
                        result['cz'] = float(bits[i + 1])
                    if bits[i] == "+/-":
                        result['czErr'] = float(bits[i + 1])
                if result['RVal'] > bestRVal:
                    result['z'] = result['cz'] / 3e5
                    result['zErr'] = (result['czErr'] / result['cz']) * \
                                     (result['cz'] / 3e5)
                    bestResult = result
                    bestRVal = result['RVal']
        return bestResult

    def loadTemplates(self, tempDir, templateFileNamesList):
        """ Takes in a list of SDSS template file names. Appends the Tremonti
        starburst template which we handle differently.

        Returns a list containing astSED.SED objects

        """

        print "Loading templates ..."

        templatesList = []
        for t in templateFileNamesList:
            # This loads in an SDSS spectral template, and feeds it into a SED
            # object
            timg = pyfits.open(tempDir + os.path.sep + t)
            th = timg[0].header
            tc0 = th['COEFF0']
            tc1 = th['COEFF1']
            tpixRange = numpy.arange(timg[0].data.shape[1])
            twavelengthRange = 10.0 ** (tc0 + tc1 * tpixRange)

            tempSED = astSED.SED(wavelength=twavelengthRange,
                                 flux=timg[0].data[0])
            templatesList.append(tempSED)

        # Tremonti star burst
        #s=astSED.SED()
        #s.loadFromFile(tremontiFileName)
        #templatesList.append(s)

        return templatesList

    def log(self):
        # Save results for current spectrum
        self.savePNG()

        self.outFile = file(outDir + os.path.sep +\
            self.objectSpecFileNames[self.currentSpecFileIndex].rstrip('.fits')
            + ".results", "a")
        #if self.qualityRadioVar.get() != 0:
        self.outFile.write("%s\t%.5f\t%.5f\t%d\t%s\n" % \
                            (self.fibernumberScale.get(),
                            self.redshiftScaleVar.get(),
                            float(self.redshiftErrorEntryVar.get()),
                            self.qualityRadioVar.get(),
                            self.commentsEntryVar.get()))
#        else:
#            self.outFile.write("%s\t%s\t%s\t%d\t%s\n" % \
#                               (self.objectSpecFileNames[
#                                   self.currentSpecFileIndex],
#                                "None", "None", self.qualityRadioVar.get(),
#                                self.commentsEntryVar.get()))

        self.outFile.close()
        print('Logged %s' % self.fibernumberScale.get())
#        # Move on to next spectrum
#        if self.currentSpecFileIndex < len(self.objectSpecFileNames) - 1:
#            self.currentSpecFileIndex += 1
#
#            print "Checking spectrum %d/%d ..." % \
#           (self.currentSpecFileIndex + 1, len(self.objectSpecFileNames))

#            obj = self.loadObjectSpectrum(self.objectSpecFileNames[
#                self.currentSpecFileIndex])
#            self.objectSpecFileName = self.objectSpecFileNames[
#                self.currentSpecFileIndex]  # for plot titles etc.
#            self.objSED = obj['object']
#            self.unsmoothedObjFlux = self.objSED.flux[:]
#            self.skySED = obj['sky']

#            x = self.smoothScale.get()
#            self.objSED.flux = self.unsmoothedObjFlux[:]
#            self.objSED.smooth(x)

#            self.resetFeatures()

#            self.qualityRadioVar.set(0)

#            self.commentsEntryVar.set("")
#            self.outPathEntryVar.set(self.outDir + os.path.sep +
#                                     self.objectSpecFileName.replace(".fits",
#                                         ".png"))
#            self.redshiftErrorEntryVar.set(0.001)

#            self.updatePlot(self.objSED,
#                            self.templates[
#                                self.templateRadioVar.get()],
#                            self.skySED,
#                            self.redshiftScaleVar.get(),
#                            tempLabel=os.path.split(templateLabels[
#                                self.templateRadioVar.get()])[-1],
#                                redrawSky=True,
#                                redrawFeatures=True,
#                                plotFeatures=self.plotFeatures)
#        else:
#            print "Finished checking all spectra!"
#            self.buttonFrame.quit()

    def loadIFUSpectra(self, objectFileName):
        """ Loads in an object spectrum - this has to be in DEEP2 pipeline
        spec1d format (i.e. fits tables)
        Object spectrum is smoothed by boxcar of size smoothPix.

        Returns a dictionary containing object and sky astSED.SED objects
        {'object', 'sky'}

        """

        print "Loading IFU spectrum ..."

        self.objectFileName = objectFileName
        oimg = pyfits.open(objectFileName)

        # Load the IFU data -- Row-stacked spectra
        odata = oimg[1].data
        odata_dim = odata.shape
        wcs = astWCS.WCS(objectFileName, extensionName=1)
        owavelengthStartEnd = wcs.getImageMinMaxWCSCoords()[0:2]
        fiberNumber = wcs.getImageMinMaxWCSCoords()[2:4]
        owavelengthStep = oimg[1].header['CDELT1']

        owavelengthRange = [owavelengthStartEnd[0] + i * owavelengthStep
                            for i in range(odata_dim[1])]

        # Check to make sure we got it right
        if not owavelengthRange[-1] == owavelengthStartEnd[-1]:
            print 'The ending wavelenghts do not match... Exiting'
            sys.exit(1)
        else:
            sums = [sum(odata[i,:]) for i in range(odata.shape[0])]
            #find the median value of all the fibers
            med = astStats.clippedMedianStdev(sums)
            med = med['clippedMedian']

            skyfibers = [i for i in range(odata.shape[0])\
                    if sum(odata[i,:]) <= med]
            skydata = odata.take(skyfibers, axis=0)

            oskyflux = [numpy.average(skydata[:,i])\
                    for i in range(skydata.shape[1])]

        RSS = []
        for i in range(int(fiberNumber[1])):
            oflux = odata[i] - oskyflux
            #oflux = odata[i]

            # Mask out extreme values in spectrum
            # Just because edges dodgy in efosc
            med = numpy.median(oflux)
            oflux[numpy.greater(abs(oflux), 10.0*med)] = 0.0001

            objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)

            #  make it > 0 everywhere
            #objSED.flux = objSED.flux - objSED.flux.min()
            #objSED.flux = objSED.flux / objSED.flux.max()

            skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)

            RSS.append({'object': objSED, 'sky': skySED})
        return RSS
        #return {'object': objSED, 'sky': skySED}

    def loadObjectSpectrum(self, objectFileName):
        """ Loads in an object spectrum - this has to be in DEEP2 pipeline
        spec1d format (i.e. fits tables)
        Object spectrum is smoothed by boxcar of size smoothPix.

        Returns a dictionary containing object and sky astSED.SED objects
        {'object', 'sky'}

        """

        print "Loading object spectrum ..."

        self.objectFileName = objectFileName
        oimg = pyfits.open(objectFileName)

        # If DEEP2 format, concatenate red, blue spectra
        # Otherwise, assume in efosc2 reducer format
        try:
            rwavelengthRange = oimg['HORNE-R'].data.field('LAMBDA')[0]
            rflux = oimg['HORNE-R'].data.field('SPEC')[0]
            rskyflux = oimg['HORNE-R'].data.field('SKYSPEC')[0]
            bwavelengthRange = oimg['HORNE-B'].data.field('LAMBDA')[0]
            bflux = oimg['HORNE-B'].data.field('SPEC')[0]
            bskyflux = oimg['HORNE-B'].data.field('SKYSPEC')[0]
            owavelengthRange = numpy.array(
                bwavelengthRange.tolist() + rwavelengthRange.tolist())
            oflux = numpy.array(bflux.tolist() + rflux.tolist())
            oskyflux = numpy.array(bskyflux.tolist() + rskyflux.tolist())
        except:
            # owavelengthRange = oimg['1D_SPECTRUM'].data.field('LAMBDA')
            # oflux = oimg['1D_SPECTRUM'].data.field('SPEC')
            # columnNames = oimg['1D_SPECTRUM'].columns.names
            # if 'SKYSPEC' in columnNames:
                # oskyflux = oimg['1D_SPECTRUM'].data.field('SKYSPEC')
            #else:
                # oskyflux = None

            owavelengthRange = oimg[1].data.field('LAMBDA')
            oflux = oimg[1].data.field('SPEC')
            columnNames = oimg[1].columns.names
            if 'SKYSPEC' in columnNames:
                oskyflux = oimg[1].data.field('SKYSPEC')
            else:
                oskyflux = None

        # Mask out extreme values in spectrum
        # Just because edges dodgy in efosc
        #med=numpy.median(oflux)
        #oflux[numpy.greater(abs(oflux), 10.0*med)]=0.0

        objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)
        #  make it > 0 everywhere
        objSED.flux = objSED.flux - objSED.flux.min()
        objSED.flux = objSED.flux / objSED.flux.max()
        if oskyflux is not None:
            skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)
        else:
            skySED = None

        return {'object': objSED, 'sky': skySED}

    def plotSkyChanged(self):
        """ Clears figure for if we're redrawing sky subplot or not. Use if we
        change norm method too.

        """
        pylab.clf()
        self.redrawPlot()

    def maskLines(self, objSED, skySED):
        # Mask prominent sky emission lines
        if skySED is None:
            fixbad = "n"
        else:
            fixbad = "y"
        if fixbad == "y":
            normSkyFlux = skySED.flux / skySED.flux.max()
            threshold = 0.25
            badPix = numpy.where(normSkyFlux > threshold)[0]
            lines = []
            for i in range(len(badPix)):
                startPixInLine = False
                for line in lines:
                    if line[0] <= badPix[i] <= line[1]:
                        startPixInLine = True
                if not startPixInLine:
                    pixCount = 1
                    if pixCount + i < len(badPix) - 1:
                        nextPix = badPix[i + pixCount]
                        prevPix = badPix[i]
                        maxReached = False
                        while nextPix < prevPix + 2:
                            if pixCount + i < len(badPix) - 1:
                                prevPix = badPix[i + pixCount]
                                pixCount += 1
                                nextPix = badPix[i + pixCount]
                            else:
                                maxReached = True
                                break
                        if not maxReached:
                            lastPix = prevPix
                        else:
                            lastPix = max(badPix)
                    else:
                        lastPix = max(badPix)
                    #append the lines and add a little padding
                    lines.append([badPix[i]-5, lastPix+5])

            for line in lines:
                # Do a simple linear fit to the end points
                y = objSED.flux[line[1]] - objSED.flux[line[0]]
                x = skySED.wavelength[line[1]] - skySED.wavelength[line[0]]
                slope = y/x
                intercept = objSED.flux[line[0]] - slope *\
                    skySED.wavelength[line[0]]
#                print skySED.wavelength[line[0]], skySED.wavelength[line[1]]

                for i in range(line[0], line[1]):
                    objSED.flux[i] = slope * skySED.wavelength[i] + intercept

            return objSED

    def updatePlot(self, objSED, tempSED, skySED, redshift, tempLabel=None,
                   redrawSky=True, redrawFeatures=False, plotFeatures=[]):
        """ Updates the pylab plot of the object spectrum with template
        overlaid.

        """

        xMin = float(self.minWavelengthEntry.get())
        xMax = float(self.maxWavelengthEntry.get())

        tempSED.redshift(redshift)
        tempSED.flux = tempSED.flux

        # Mask out the spectrum that we can't see
        plusMask = numpy.greater(tempSED.wavelength, xMin)
        minusMask = numpy.less(tempSED.wavelength, xMax)
        mask = numpy.logical_and(plusMask, minusMask)
        tempSED.flux = tempSED.flux - tempSED.flux.min()
        tempSED.flux = tempSED.flux / tempSED.flux[mask].max()

        # In case we don't want to see redward of 10000 Angstroms
        plusMask = numpy.greater(objSED.wavelength, xMin)
        minusMask = numpy.less(objSED.wavelength, xMax)
        mask = numpy.logical_and(plusMask, minusMask)

        objSED = self.maskLines(objSED, skySED)

        # Normalize
        objSED.flux = objSED.flux - objSED.flux[mask].min()
        objSED.flux = objSED.flux / objSED.flux[mask].max()

        # Norm based on matching flux as closely as possible between template
        # and object spectrum
        # Ignore XX% of each end of spectrum as edges can be weird
        if self.altNormCheckVar.get() == 1:
            ignoreAngstroms = (objSED.wavelength.max() -
                               objSED.wavelength.min()) * 0.25
            dw = 100
            binEdges = numpy.arange(objSED.wavelength.min() + ignoreAngstroms,
                                    objSED.wavelength.max() - ignoreAngstroms,
                                    dw)
            passbands = []
            for b in binEdges:
                passbands.append(astSED.TopHatPassband(b, b + dw))
            objSEDDict = objSED.getSEDDict(passbands)
            tempSEDDict = tempSED.getSEDDict(passbands)
            norm0 = 1.0
            norm, success = optimize.leastsq(fitSEDNormErrFunc, norm0,
                                             args=(tempSEDDict['flux'],
                                                   objSEDDict['flux']))
            objSED.flux = objSED.flux / norm

            #if skySED != None and self.plotSkyCheckVar.get() == 1:
            #pylab.subplot(211)
        pylab.cla()
        pylab.title(self.objectSpecFileName+' Fiber No. ' +
                str(self.fibernumberScale.get()))
        pylab.plot(objSED.wavelength, objSED.flux, 'k-')
        pylab.plot(tempSED.wavelength, tempSED.flux, c='#A60628',
                label=tempLabel)
        pylab.text(0.05, 0.92, "z = %.5f $\pm$ %.5f (Q = %s)" %
                               (tempSED.z,
                                float(self.redshiftErrorEntryVar.get()),
                                self.qualityRadioVar.get()),
                   ha='left',
                   va='top',
                   transform=pylab.gca().transAxes, size=12, color='#A60628')

        pylab.ylim(0, 1.2)

        #if skySED != None and self.plotSkyCheckVar.get() == True:
        #pylab.gca().set_xticklabels([])
        #else:
        pylab.xlabel("Wavelength (Angstroms)")

        # Plots the spectral features in turn
        #plotFeatures=["H", "K", "[OII]"]
        if redrawFeatures:
            #ylim=pylab.gca().get_ylim() # Need this to automatically draw
            #correct length -- for features
            for on, item in zip(self.plotFeatures, spectralFeaturesCatalogue):
                if on == 1:
                    featureLabel = item[0]
                    # Greek letters? eta will cause a problem here!
                    featureLabel = featureLabel.replace("alpha",
                                                        "$\\alpha$")
                    featureLabel = featureLabel.replace("beta",
                                                        "$\\beta$")
                    featureLabel = featureLabel.replace("gamma",
                                                        "$\gamma$")
                    featureLabel = featureLabel.replace("delta",
                                                        "$\delta$")
                    featureLabel = featureLabel.replace("epsilon",
                                                        "$\\epsilon$")
                    featureLabel = featureLabel.replace("zeta",
                                                        "$\zeta$")
                    featureLabel = featureLabel.replace("theta",
                                                        "$\\theta$")
                    for i in range(1, len(item)):
                        featureLambda = (1.0 + float(redshift)) * item[i]
                        pylab.plot((featureLambda, featureLambda), (0, 1.0),
                                   '--', c='#467821')
                        pylab.text(featureLambda, 1.05, featureLabel,
                                   ha='center', va='top', size=10,
                                   rotation='vertical')

        if redrawSky and self.plotSkyCheckVar.get() == 1:
            if skySED is not None:
                pylab.plot(skySED.wavelength,
                           skySED.flux / skySED.flux.max() * 0.3, c='#348ABD',
                           label='Sky')
            # Main telluric absorption features
            c = patches.Rectangle((6860, 0), (6930 - 6860), 1.2, fill=True,
                                  edgecolor=(0.8, 0.8, 0.8),
                                  facecolor=(0.8, 0.8, 0.8), linewidth=1)
            pylab.gca().add_patch(c)
            c = patches.Rectangle((7590, 0), (7710 - 7590), 1.2, fill=True,
                                  edgecolor=(0.8, 0.8, 0.8),
                                  facecolor=(0.8, 0.8, 0.8), linewidth=1)
            pylab.gca().add_patch(c)

        # Finish drawing the object spectrum plot
        pylab.ylim(0, 1.2)
        pylab.xlim(xMin, xMax)
        pylab.ylabel("Relative Flux")
        pylab.xlabel("Wavelength (Angstroms)")
        pylab.legend(loc="upper right")

    def redrawPlot(self):
        self.updatePlot(self.obj[246-self.fibernumberScale.get()]['object'],
                self.templates[self.templateRadioVar.get()],
                self.obj[246-self.fibernumberScale.get()]['sky'],
                self.redshiftScaleVar.get(),
                tempLabel=os.path.split(templateLabels[
                    self.templateRadioVar.get()])[-1],
                redrawSky=True,
                redrawFeatures=True,
                plotFeatures=self.plotFeatures)
#        self.updatePlot(self.objSED,
#                        self.templates[self.templateRadioVar.get()],
#                        self.skySED,
#                        self.redshiftScaleVar.get(),
#                        tempLabel=os.path.split(templateLabels[
#                            self.templateRadioVar.get()])[-1],
#                        redrawFeatures=True,
#                        plotFeatures=self.plotFeatures)


def fitSEDNorm(p, modelFluxes):
    """ Pair of helper functions for fitting SED normalisation
    p0 is list, [0] = normalisation

    """
    result = p * modelFluxes
    return result


def fitSEDNormErrFunc(p, modelFluxes, observedFluxes):
    x = fitSEDNorm(p, modelFluxes) - observedFluxes
    chiSq = numpy.sum(x ** 2)   # not really chi sq, duh
    return chiSq

#-----------------------------------------------------------------------------
# Main ...
if __name__ == "__main__":

    if len(sys.argv) < 3:

        print "Run: visualTemplateRedshift5.py <spec1d object spectra .fits", \
            "[wildcards allowed]> ... <outputDir>"

    else:

        objectSpecFileNames = sys.argv[1:-1]
        outDir = sys.argv[-1]

        print "File to be used: "
        print objectSpecFileNames

        root = Tkinter.Tk()
        root.title("Visual Template Redshift 5.0 now with 100% more IFU")
        app = App(root, objectSpecFileNames, outDir)
        root.mainloop()
