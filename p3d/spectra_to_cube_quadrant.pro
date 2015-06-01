;#############################################################################
;
; Based initially off code written by Sarah Brough, Australian Astronomical Observatory
;
; Last updated by Jimmy
; E-mail: jimmy@physics.tamu.edu
; 
; Updated versions of the software are available from my web page
; http://galaxies.physics.tamu.edu/jimmy/
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;
; NAME:
;   SPECTRA_TO_CUBE_QUADRANT
;
; PURPOSE:
;   This code creates a data cube for each quadrant containing the spectra.
;   Variance data and spectra with skylines are added as seperate extensions.
;	During the proces it scales for differences in fiber transmission and
;	subtracts the sky
;
;
; CALLING SEQUENCE:
;   spectra_to_cube_quadrant,'quadrant number','observation','pointing'
;	eg spectra_to_cube_quadrant,'1','a','pt1'
;
; INPUT PARAMETERS:
;   Quadrant Number: A number 1-4 identifying the Quadrant, starting with 1 in the 
;		upper right corner and going counter-clockwise.
;   Observation: Galaxies have multiple pointings, this identifies which
;		pointing the data is from, identified chronologically by letter
;   Pointing: To identify which pointing we're reducing, different pointings
;		need different calibration files, and may be rotated by 180 degrees
;		which would be identified as pt180, pt181, etc.
;
; ENVIRONMENTAL VARIABLES:
;	If called by a bash script, the following variables must be defined in the bash
;	script that called this program.
;
;	infile1: The directory containing the spectra files produced by the pipeline.
;	infile2: The calibration files directory.
;	fiber_mask: Filename of the list of bad fibers.
;	'skyfiber'+Quad+obs+'start': The first fiber that defines the sky fiber region
;	'skyfiber'+Quad+obs+'end': The last fiber that defines the sky fiber region
;
; OUTPUT:
;   One data cube with 3 extentions.  Skylines are kept in one extention to be
;		used later when normalizing transmission across the whole cube.
;
; NOTES:
;	If run directly from IDL, edit everything within an 'if (testing ne 1)'
;		statement to have the proper directories.
;	May need to edit wave_pix1 & wave_pix2 depending on instrument and pipeline.
;		These parameters determine the pixels used for background noise
;		determination and should be taken from a quiet part of the spectrum near
;		the beginning. Currently set for the VIMOS pipeline.
;
;--------------------------------
;
; LOGICAL PROGRESSION OF IDL CODE:
;	1.Read in the data, and reference files.
;	2.Set Bad fibers to zero
;	3.Determine Variance
;	4.Determine Sky Fibers (Sky fibers are also used for noise determination)
;	5.Determine Noise
;	6.Calculate and Subtract Noise
;	7.Determine Gaussian Fit to Sky Fibers
;	8.Normalize Transmission
;	9.Remove Sky
;	10.Create Cubes
;	11.Write Cubes
;
;--------------------------------
;
; REQUIRED ROUTINES:
;       IDL Astronomy Users Library: http://idlastro.gsfc.nasa.gov/
;		Heirarch.pro http://www.astro.uu.se/~piskunov/RESEARCH/REDUCE/reduce.html
;		Hitme.pro Written by Rob Sharp Australian Astronomical Observatory
;
; MODIFICATION HISTORY:
;   V1.0 -- Created by Jimmy, 2011
;
;----------------------------------------------------------------------------

pro spectra_to_cube_quadrant,quad,obs,pointing

testing=0 ;Set to 0 if you want to be in "testing" mode, where more output is displayed, and files are split up, so they can be more easily examined, also paths are then hard coded.
testing_string=getenv('not_testing') ;if called from a bash script, these two lines pull it out of testing mode
testing=FIX(testing_string)

;Defined way ahead of time because these 3 lines are most edited when checking sky fibers
if (testing ne 1) then begin
    sky_fiber_1=242
    sky_fiber_2=278
    gal_name='2039'
endif

;READ IN THE DATA AND REFERENCE FILES
;Assign the filename to the variable file, then read in using mrdfits from NASA

if (testing) then begin
;environmental variables are set with the export command in bash.
    file=getenv('infile1')+'ifu_science_reduced'+Quad+''+obs+'_'+pointing+'.fits'
    ifu_table_file=getenv('infile2')+'ifutableHR'+Quad+'.fits'
endif

if (testing ne 1) then begin
        file='/Users/jimmy/Astro/reduced/'+gal_name+'pro/ifu_science_reduced'+Quad+''+obs+'_'+pointing+'.fits'
        ifu_table_file='/Users/jimmy/Astro/reduced/cal/ifutableHR'+Quad+'.fits'
endif

;input: (filename,extention number,header), to read in VIPGI produced data, change extention from 0 to 10
;result: (combined image array with (x= wavelength, y = fibre number)
im=mrdfits(file, 0, header0) ;read in first extension which contains the raw image data (spectra by fiber number)
ifu_table=mrdfits(ifu_table_file,1,header1) ;read in ifu binary table with information on where each fiber is located in the final images/telescope design

;we'll need to know how big our arrays are, so we create a variable imsize to store those dimensions.
imsize=size(im)
number_of_wave_pixels=imsize[1]
number_of_fibers=imsize[2]

;These tell us about the scales of our data, so we're going to want to read them in and put them into the final cubes
crval=fxpar(header0,'CRVAL1')
crpix=fxpar(header0,'CRPIX1')
cdelt=fxpar(header0,'CDELT1')

if (crpix le 1) then crpix=0.0 else crpix=crpix-1 ;set our reference pixel to 0 because idl will count from 0 in all our arrays

if (testing ne 1) then begin
    print,'crval:',crval
    print,'crpix:',crpix
    print,'cdelt:',cdelt
endif



;SET BAD FIBERS TO ZERO
;Done by importing a single column list of bad fibers (with a 0 at the end to keep the loop from freaking out)

bad=fltarr(number_of_wave_pixels,number_of_fibers) ;fibers listed as bad
zero=fltarr(number_of_wave_pixels,number_of_fibers) ;fibers with zero flux, both natural and the ones listed as bad
bright=fltarr(number_of_wave_pixels,number_of_fibers) ;fibers that are super bright, probably bad for another reason
fiber_ident=fltarr(number_of_wave_pixels,number_of_fibers) ;used in the diagnostic cube to identifiy which fiber is where in the field of view
bad_fibers = findgen(number_of_fibers) ;stores the column list that gets read in

if (testing ne 1) then begin
    readcol,'/Users/jimmy/Astro/reduced/'+gal_name+pointing+'sof/bad_fibers_'+Quad+'.txt',bad_fibers, FORMAT='F'
endif

if (testing) then begin
    readcol,getenv('fiber_mask'),bad_fibers, FORMAT = 'F'
endif

;Bad and dead fiber setting and fiber analysis stuff.
h=0 ;because the length of our bad fibers list is not 400 fibers long, we use h as a tracker
for i=0,number_of_fibers-1 do begin
    fiber_ident[*,i]=i+1 ;useful for identifying fiber number if something wrong with one of them
    if ((bad_fibers[h]-1) eq i) then begin
        im[*,i]=0.0 ;if the fiber is bad, set it to zero so it doesn't mess things up.
        bad[*,i]=i+1 ;set the bad fiber to it's number in the diagnostic cube so it can be identified if needed
        h=h+1
    endif
    if (max(im[*,i]) le 0.0) then begin
        zero[*,i]=i+1 ;Creates a list of what fibers are zero, useful for checking what's dead and what's been set to zero.
    endif
endfor


;if pointing eq '1153' then begin 
;	;Smoothing used for testing, do not keep this around forever.
;	new_im = fltarr(number_of_wave_pixels,number_of_fibers)
;	for a=1, 10 do begin
;		for b=1, 20 do begin
;			print,(40*(a-1))+b-1,'+',(40*(a-1))+b,'+',(40*a)-b,'+',(40*a)-b-1
;			new_im[*,(40*(a-1))+b-1]=(im[*,(40*(a-1))+b-1]+im[*,(40*(a-1))+b]+im[*,(40*a)-b]+im[*,(40*a)-b-1])/4
;			new_im[*,(40*(a-1))+b]=(im[*,(40*(a-1))+b-1]+im[*,(40*(a-1))+b]+im[*,(40*a)-b]+im[*,(40*a)-b-1])/4
;			new_im[*,(40*a)-b]=(im[*,(40*(a-1))+b-1]+im[*,(40*(a-1))+b]+im[*,(40*a)-b]+im[*,(40*a)-b-1])/4
;			new_im[*,(40*a)-b-1]=(im[*,(40*(a-1))+b-1]+im[*,(40*(a-1))+b]+im[*,(40*a)-b]+im[*,(40*a)-b-1])/4
;			b=b+1
;		endfor
;		;print,'SMOOTHING IS DONE'
;	endfor
;	
;	for c=0, 399 do begin
;		im[*,c]=new_im[*,c]
;	endfor
;endif

if (testing ne 1) then begin
    ;This is used to test what the image files look like when you black out the fibers from the bad fiber mask.
    image_out='/Users/jimmy/Downloads/fiber_lines.fits'
    mwrfits,im,image_out,create=1
    head=headfits(image_out) ;reads header from newly created file so that does not create header with wrong structure
    ;add important parameters to the header
    sxaddpar,head,'CRVAL3',crval
    sxaddpar,head,'CRPIX3',crpix
    sxaddpar,head,'CDELT3',cdelt
    mwrfits,im,image_out,head,create=1 ;write over former cube with new header
endif




;SKY FIBER SELECTION
;Select a region of the sky that will be considered the sky fibers, the whole range will be used and averaged together.
;Sky fibers are pulled from the bash script or defined above if testing.

if (testing) then begin
    sky_fiber_1=FIX(getenv('skyfiber'+'_'+Quad+obs+'_'+pointing+'_'+'start'))
    sky_fiber_2=FIX(getenv('skyfiber'+'_'+Quad+obs+'_'+pointing+'_'+'end'))
endif

;idl counts from zero whereas the VIMOS pipeline counts from 1.
sky_fiber_1 = sky_fiber_1 - 1
sky_fiber_2 = sky_fiber_2 - 1




;DETERMINE THE VARIANCE
;Variance is calculated using the first portion of the sky fibers, which are assumed to have no image data, and the first section should have no absorbtion features
;Noise model is given by: error = sqrt(electrons)/gain = sqrt(count*gain)/gain = sqrt(count/gain). Then add in some readout noise, estimated from the blank bits of the array

gain=hierarch(headfits(file), 'HIERARCH ESO DET OUT1 GAIN') ;use the hierarch.pro file I found on the internet to pull out the gain from the hierarchical header structure
if (testing ne 1) then begin
    print, 'Gain read in:',gain ; Electorns per ADU (actually inverse gain)
endif

fiber=findgen(number_of_fibers) ;Used to when plotting to be the fiber axis.

;In order to find the variance, we pick a portion of the spectrum that should just be noise
wave_pix1=1100 ;beginning of the spectrum that is just noise, changes if cdelt changes.
wave_pix2=1250 ;set our window to be 150 wavelength pixels wide

if (testing ne 1) then begin
    ;This is our check of what's actually noise
    window,4,xsize=1200,ysize=700
    plot,fiber,total(im[wave_pix1:wave_pix2,*],1),yrange=[0.0,30000.0],xtitle="fibre number",ytitle="summed counts",title="Nrd Determination" ;summed counts is between assigned wavelengths
    oplot,[sky_fiber_1,sky_fiber_1],!y.crange,linestyle=2 ;dotted lines show the region we are using
    oplot,[sky_fiber_2,sky_fiber_2],!y.crange,linestyle=2
    hitme ;used to check if this is a good noise region while testing variance
endif

Nrd=robust_sigma(im[wave_pix1:wave_pix2,sky_fiber_1:sky_fiber_2])
if (testing ne 1) then begin
    print,'Estimated standard deviation:',Nrd
endif

;Create the variance array, The FLTARR function creates a floating-point vector or array of the specified dimensions.
var1=fltarr(number_of_wave_pixels,number_of_fibers)

;calculate variance
for i=0,number_of_fibers-1 do begin
    ;Sarah had the first version, the second version is actual variance, but I like the third
    ;var1[*,i] = abs(im[*,i]/gain)+ Nrd^2 ;variance is the gain adjusted image frame + the standard deviation squared
    ;var1[*,i] = Nrd^2 ; the standard deviation squared, variance
    var1[*,i] = Nrd ; the standard deviation, this is what I like, plugs into pPXF well.
endfor




;CALCULATE AND SUBTRACT THE NOISE FROM OUR SIGNAL
;Using the sky fibers, and the signal in the 5577 angstrom region, we calculate the vertical offset, which is just background noise, and subtract that away.
 
;Create a vector listing the wavelengths, then select only those wavelengths from that vector that are in the 5577 skyline region, will have final answer in pixels.
wave_pix_test=range(0,number_of_wave_pixels-1) ;Creates a vector, 4259 elements big, each element is set to it's own wavelength pixel number, or wave_pix_test[n]=n
wave_angst_test=crval+((wave_pix_test-crpix)*cdelt) ;Each element in the array is converted to it's wavelength in angstroms
test = bytarr(number_of_wave_pixels) ;The BYTARR function creates a byte vector or array.  Everything set to zero.
test or= wave_angst_test gt 5557 and wave_angst_test lt 5597 ;sets the test vector to 1 when in the skyline range, 0 otherwise
in_pix_range=where(test eq 1) ;pixel numbers that are within the skyline wavelength range
wave_pix_skyrange=wave_pix_test(in_pix_range) ;sets the range for the skyline in pixels
 
im_corr=fltarr(number_of_wave_pixels,number_of_fibers) ;float array to represent the image with offset subtracted

;calculate noise around 5577 skyline, that's the offset, found using a robust linear fit
for i=0,number_of_fibers-1 do begin
    ;This if statement stops the "CURVEFIT: Failed to converge- CHISQ increasing without bound" errors when there's nothing there to fit.
    if (max(im[in_pix_range,i]) gt 1) then begin
 
        test=robust_linefit(wave_pix_skyrange,im[in_pix_range,i],offset) ;test is a throw away variable here that holds the coefficient, COEFF = ROBUST_LINEFIT( X, Y, YFIT)
 
        im_corr[in_pix_range,i]=im[in_pix_range,i]-offset ;subtract the noise from the image array, assign it to a new array that should be zero when there's no signal.
        
        if (testing ne 1) then begin    
            ;plotting checks, the white straight line should lie along the horizontal aspects of the signal.
            wave_pix_5557 = ((5557-crval)/cdelt)+crpix ;for the 5577 skyline, convert to pixels
            wave_pix_5597 = ((5597-crval)/cdelt)+crpix ;set up a 2nd limit, 40 angstroms away.
            plot,wave_pix_skyrange,im[wave_pix_5557:wave_pix_5597,i],yrange=[-10.,4000.],xtitle="pixel",ytitle="counts" ;solid curved line, is the original signal line
            oplot,wave_pix_skyrange,offset,linestyle=3 ;straight horizontal line, this is the noise that has since been removed
            legend,'Fiber Number'+STRING(i)

            print,i ;slows down the plotting so we can actually see what's going on.
            ;hitme ;to stop things all together
        endif
     
    endif
endfor
 
if (testing ne 1) then begin    
    hitme ;to stop between noise subtraction and flux calculation.
endif
 
  
;DETERMINE FLUX NORMALIZATION
;Perform a Gaussian fit on the data with the offset removed, and then come up with a normalization factor for each fiber.
;Using the skyline because it's very prominant and easy to fit to, should also be the same in each fiber, so we can use that to normalize everything.
 
sky_flux_5577=fltarr(number_of_fibers) ;integrated flux in skyline
 
for i=0,number_of_fibers-1 do begin
    if (max(im[in_pix_range,i]) gt 1) then begin
        ;using mpfitpeak removes the floating underflow errors, but doesn't fit right.
        skyfit=gaussfit(wave_pix_skyrange,im_corr[in_pix_range,i],gauss_coeff) ;fit a guassian to the skyline with noise subtracted

        if (testing ne 1) then begin 
            ;plotting checks, the similarly colored lines should match up very closely, that is if the fit is matching the signal properly.
            plot,crval+((wave_pix_skyrange-crpix)*cdelt),im[wave_pix_5557:wave_pix_5597,i],yrange=[-10.,3500.],xtitle="pixel",ytitle="counts" ;solid curved line, is the original signal line
            oplot,crval+((wave_pix_skyrange-crpix)*cdelt),skyfit,color=175, linestyle=2 ;dotted line, gaussian fit to the continuium
            oplot,crval+((wave_pix_skyrange-crpix)*cdelt),im_corr[in_pix_range,i],color=175 ;this is the line that we were trying to fit to
            legend,'Fiber Number'+STRING(i)
        endif
 
        ;calculate integrated flux in sky line
        sky_flux_5577[i]=gauss_coeff[0]*gauss_coeff[2]*sqrt(2*!pi) ;this is what's in the skyline, used to normalize the rest of the data
        
        if (testing ne 1) then begin
            print,'fiber number',i,' Mean of 5577 skyline fit:',crval+((gauss_coeff[1]-crpix)*cdelt) ;check that the mean value stored in mean_5577 is close to the 5577.33 skyline
            ;hitme ;if needed to go 1 by 1, printing above also does a good job of slowing things down.
        endif
    endif
endfor
 
 
med_int_flux=median(sky_flux_5577) ;calculate the median of the flux in all fibers
if (testing ne 1) then begin
    print, "Median Flux in 5577 skyline: ",med_int_flux ;note to self, not a wavelength, but an intensity
endif
scale=sky_flux_5577/med_int_flux ;divide lines by median integrated flux so all are compared to median, so fiber flux in the skyline is now relative to the mid point

for i=0,number_of_fibers-1 do begin
    if (scale[i] gt 1.5) then begin
        bright[*,i]=i+1
    endif
endfor
 
if (testing ne 1) then begin
    ;checking integrated flux measurements
    plot,fiber,scale,xtitle="fibre number",ytitle="Scale factor",yrange=[0,1.5] ;this is the fiber transmission value plotted by fiber. 
    hitme ;if you need to check the ratios of the scale.
endif
 
 
 
;NORMALIZE TRANSMISSION
;Normalize the transmission for each fiber relative to the median.
 
im_trans=fltarr(number_of_wave_pixels,number_of_fibers) ;will be used to store the transmission normalized image data
var_trans=fltarr(number_of_wave_pixels,number_of_fibers) ;will be used to store the transmission normalized variance data
 
;normalize fiber transmission from data & variance
for i=0,number_of_fibers-1 do begin
    ;Properly handle the cases where the scale is zero, don't want to divide by zero.
    if (scale[i] gt 0) then begin
        im_trans[*,i]=(im[*,i]/scale[i]) ;image transmission is the uncorrected image data scaled by our normalization
        var_trans[*,i]=(var1[*,i]/(scale[i])^2) ;as is multiplicative need to adjust variance too
    endif else begin
        im_trans[*,i]=0.0 ;if the flux is zero, then the the im_trans should also be zero
        var_trans[*,i]=9.9e6^2 ;zero variance would mean good data, so set it very high.
    endelse
endfor
 
 
 
 
;DETERMINE SKY FOR SUBTRACTION
;Create a sky vector for each fiber, containing the intensities at that pixel in that region that is skyline, and subtract it from transmission scaled data.
 
sky_vect=fltarr(sky_fiber_2-sky_fiber_1+1) ;1 element larger than sky pixels, holds the sky intensity before it gets averaged for subtraction
sky=fltarr(number_of_wave_pixels) ;sky is the actual sky intensity that gets subtracted from the normalized image data.
 
for i=0,number_of_wave_pixels-1 do begin ;i is the fiber tracker
    j=0 ;j is used as a pixel tracker
    
    ;create the sky vector array for this fiber
    for k=sky_fiber_1,sky_fiber_2 do begin
        ;if the total of the intensity (in the wavelength dimension) is not zero ,1 chooses the dimension
        if (total(im_trans[i,k],1) gt 0.0) then begin            
            sky_vect[j]=im_trans[i,k] ;sky vector is the intensity of the pixels for that fiber in the skyline region
            j=j+1 ;increase j to set the next wavelength pixel to cover the whole skyline region
        endif        
    endfor
    
    ;if fiber isn't 0 then calculate the mean and store the average as our actual sky
    if (total(sky_vect[*],1) gt 0.0) then begin
        meanclip, sky_vect, output_mean,clipsig=5 ;Compute the sigma clipped mean at 5 sigma on the sky vector and stores it as output_mean
        sky[i]=output_mean
    endif
    
    ;reset the vector being re-used.
    for n=0,j-1 do begin
        sky_vect[n]=0.0
    endfor
    
endfor
 
wavelength=findgen(number_of_wave_pixels) ;array of the wavelength pixels
wavelength=(wavelength*cdelt)+crval ;convert to angstroms

if (testing ne 1) then begin
    plot,wavelength,sky,xtitle="Wavelength",ytitle="Counts",title="Sky signal to be subtracted",yrange=[-10.,200.] ;plot of the sky signal that we will be later subtracting
    hitme ;If we need to check what the sky signal looks like.
endif
 
im_trans_sky=fltarr(number_of_wave_pixels,number_of_fibers) ;image data, transmission and skyline corrected
 
for i=0,number_of_fibers-1 do begin 
    im_trans_sky[*,i]=im_trans[*,i]-sky[*] ;subtract the skyline from our transmission scaled image file

    if (testing ne 1) then begin    
        ;Most important plot check, this should tell us whether or not our skyline fit matches the actual sky line
        plot, wavelength,im_trans[*,i],xtitle="Wavelength",ytitle="Counts",title="Signal Subtraction Check",xrange=[5557,5597] ;Original signal
        oplot, wavelength,im_trans_sky[*,i],color=red ;New Signal
        oplot, wavelength,sky,color=blue ;subtracted signal
    
        print, i ;using this print statement slows everything down.
        ;hitme ;use hitme to really slow everything down
    endif
    
    ;this resets our final images to zero, because we subtracted above, and might have a negative flux
    if (total(im_trans[*,i],1) le 0.0) then begin     
        im_trans_sky[*,i]=0.0
    endif
endfor

sky_fibers=fltarr(number_of_wave_pixels,number_of_fibers)
for i=sky_fiber_1,sky_fiber_2 do begin
    sky_fibers[*,i]=i+1
endfor
 
 
 

;CREATE DATA CUBES 
;Import information from binary table, and use it to position the fibers in the field of view style images

;Pull in the x and y coordinates
x1=min(ifu_table.l,max=x2)
y1=min(ifu_table.m,max=y2)

;Set the actual size of the cube
Nx=max(ifu_table.l)
Ny=max(ifu_table.m)
if (testing ne 1) then begin
    print,'Size of cube : ',Nx,Ny,Nx*Ny,number_of_wave_pixels
endif

;create the cubes to be filled
cube=fltarr(Nx,Ny,number_of_wave_pixels) ;final corrected cube cube
cube_nosky=fltarr(Nx,Ny,number_of_wave_pixels) ;corrected cube
var=fltarr(Nx,Ny,number_of_wave_pixels) ;variance cube
diag=fltarr(Nx,Ny,5) ;cube of diagnostic data

;populate the cube
for i=0,number_of_fibers-1 do begin
    n1=ifu_table[i].l-1 ;x
    n2=ifu_table[i].m-1 ;y
    cube[n1,n2,*]=im_trans_sky[*,i] ;cube is the final cube with all corrections
    cube_nosky[n1,n2,*]=im_trans[*,i] ;cube_nosky is the cube without the sky subtraction
    var[n1,n2,*] = var_trans[*,i] ;var is the variance cube
    ;if (testing) then begin
        diag[n1,n2,0]=fiber_ident[0,i]
        diag[n1,n2,1]=sky_fibers[0,i]
        diag[n1,n2,2]=bad[0,i]
        diag[n1,n2,3]=zero[0,i]
        diag[n1,n2,4]=bright[0,i]
    ;endif
endfor




;WRITING INDIVIDUAL QUADRANT OUT AS CUBE

if (testing) then begin
    fileout=getenv('infile1')+'ifu_science_reduced_'+Quad+obs+'_idl_'+pointing+'.fits'
    diagout=getenv('infile1')+'diagnostic_'+Quad+obs+'_'+pointing+'.fits'
endif

if (testing ne 1) then begin
    fileout='/Users/jimmy/Downloads/cube.fits'
    fileout2='/Users/jimmy/Downloads/cube_with_sky.fits'
    fileout3='/Users/jimmy/Downloads/cube_var.fits'
    diagout='/Users/jimmy/Downloads/diagnostic.fits'
endif

print,'Writing test data cube to file : '+fileout

mwrfits,cube,fileout,create=1

head=headfits(fileout) ;reads header from newly created file so that does not create header with wrong structure
sxaddpar,head,'CRVAL3',crval ;add important parameters to the header
sxaddpar,head,'CRPIX3',crpix
sxaddpar,head,'CDELT3',cdelt

if (testing) then begin
    mwrfits,cube,fileout,head,create=1 ;write over former cube with new header, create 1 creates new file even if old one exists
    mwrfits,cube_nosky,fileout,head,create=0 ;add as new extension cube without sky_subtraction array
    mwrfits,var,fileout,head,create=0 ;add as new extension variance array
    mwrfits,diag,diagout,head,create=1 ;write our diagnostic quadrant out
endif

if (testing ne 1) then begin
    mwrfits,cube,fileout,head,create=1 ;write over former cube with new header, create 1 creates new file even if old one exists
    mwrfits,cube_nosky,fileout2,head,create=0 ;add as new extension cube without sky_subtraction array
    mwrfits,var,fileout3,head,create=0 ;add as new extension variance array
    mwrfits,diag,diagout,head,create=1 ;write our diagnostic quadrant out
endif



end