#!/bin/bash

################################################################################
# Declare variables
# These can be overwritten by a configfile given as first parameter to
# the reduction script
################################################################################

# Directory names
BIAS="zero"
FLAT="flat"
COMP="comp"
SCI="object"
EVE="eve"
MORN="morn"

# Tool paths
FT=$CUREBIN
FT2=$CUREBIN

# Data directories
#LINES=../lines.may13 
LINES=../newlines.may13 
DATADIR=./
SUBDIRS=""

# Namebase for file to be reduced
NAMEBASE="bcs"
#NAMEBASE="vp"

# If empty these values are taken from FFILE_?
BIASSEC_L=""
BIASSEC_S=""
TRIMSEC_L=""
TRIMSEC_S=""

# Configure steps to be run
MAKE_BIAS=1
MAKE_FLAT=1
MAKE_DEFORM=1
MAKE_SKYSUB=1

FIX_H2=0
FIX_H3=1

BIASLOGFILE=bias.log
FLATLOGFILE=flat.log
DEFORM1LOGFILE=deform1.log
DEFORM2LOGFILE=deform2.log
SKYLOGFILE=sky.log
ERRORLOGFILE=error.log
SCRIPTLOGFILE=script.log

LOGPATH=$DATADIR

DEBUG=false
LOGSCRIPT=false

# Additional options for generating the distortion & fibermodel on the flats
DEFORM1_OPTS=""

# Additional options for adjusting the distortion & fibermodel for each science frame
#    (add -S 1000 for HET data to force a high S/N limit on adjusting the flats to each science frame
#     effectively turning off the generation of adjusted Distortion and FiberModel objects.)
DEFORM2_OPTS=""

# Additional options for sky subtraction (add -J for HET data: force using the model flat)
SUBSKY_OPTS=""



################################################################################
# End of user definable data
################################################################################

# Check that the necessary environment variables are set
if [ x$CUREBIN == "x" ]; then
    echo "Please set CUREBIN environment variable. e.g."
    echo "bash : export CUREBIN=~/cure/bin
    echo "csh  : setenv CUREBIN "~/cure/bin"
	exit 1
fi

# Check that subtractfits is in place (we will assume that the rest
#  of the fitstools is as well if subtractfits is.
if [ ! -x $CUREBIN/subtractfits ]; then
    echo "Unable to find $FT/subtractfits, please check that CUREBIN is"
    echo "  set correctly."
    exit 1
fi

# Check that deformer3 is in place (we will assume that the rest
#  of the cure binaries is as well if deformer3 is.
if [ ! -x $CUREBIN/deformer3 ]; then
	echo "Unable to find $FT/deformer3, please check that CUREBIN is"
	echo "  set correctly."
	exit 1
fi


# Source the configfile
if [ $# -gt 0 ]; then
    . $1
fi

PWD=`pwd`

H2=""
H3=""

if [ $FIX_H2 -eq 1 ]; then
	H2="-q"
fi

if [ $FIX_H3 -eq 1 ]; then
	H3="-Q"
fi
	
BIASLOG=${PWD}/${LOGPATH}/$BIASLOGFILE
FLATLOG=${PWD}/${LOGPATH}/$FLATLOGFILE
DEFORM1LOG=${PWD}/${LOGPATH}/$DEFORM1LOGFILE
DEFORM2LOG=${PWD}/${LOGPATH}/$DEFORM2LOGFILE
SKYLOG=${PWD}/${LOGPATH}/$SKYLOGFILE
ERRORLOG=${PWD}/${LOGPATH}/$ERRORLOGFILE
SCRIPTLOG=${PWD}/${LOGPATH}/$SCRIPTLOGFILE

# figure out stuff about the data
#FFILE=`ls -1 $BIAS/${NAMEBASE}*.fits | head -1`
FFILE_L=`ls -1 $DATADIR/$BIAS/${NAMEBASE}*.fits | head -1`
FFILE_S=`ls -1 $DATADIR/$BIAS/${NAMEBASE}*.fits | head -1`

if [ "XBIASSEC_L" == "X" ]; then
    echo "Getting BIASSEC_L from $FFILE_L"
    BIASSEC_L=`$FT/headfits -r -k BIASSEC $FFILE_L | sed 's/[]\[]//g'`
fi
if [ "XBIASSEC_S" == "X" ]; then
    echo "Getting BIASSEC_S from $FFILE_S"
    BIASSEC_S=`$FT/headfits -r -k BIASSEC $FFILE_S | sed 's/[]\[]//g'`
fi
if [ "XTRIMSEC_L" == "X" ]; then
    echo "Getting TRIMSEC_L from $FFILE_L"
    TRIMSEC_L=`$FT/headfits -r -k DATASEC $FFILE_L | sed 's/[]\[]//g'`
fi
if [ "XTRIMSEC_S" == "X" ]; then
    echo "Getting TRIMSEC_S from $FFILE_S"
    TRIMSEC_S=`$FT/headfits -r -k DATASEC $FFILE_S | sed 's/[]\[]//g'`
fi

# Clean the logfiles
if [ $MAKE_BIAS -eq 1 ]; then
    cat /dev/null > $BIASLOG
fi
if [ $MAKE_FLAT -eq 1 ]; then
    cat /dev/null > $FLATLOG
fi
if [ $MAKE_DEFORM -eq 1 ]; then
    cat /dev/null > $DEFORM1LOG
    cat /dev/null > $DEFORM2LOG
fi
if [ $MAKE_SKYSUB -eq 1 ]; then
    cat /dev/null > $SKYLOG
fi

cat /dev/null > $ERRORLOG

function isbinned
{
    size=`ls -lL $1 2> /dev/null | tr -s ' ' | cut -d ' ' -f 5`
    if [ $size -eq 8527680 ]; then 
        return 0
    else
        return 1
    fi
        
}

function tolog
{
    if [ "X$LOGSCRIPT" != "Xfalse" ]; then
        echo Running: $@ >> $SCRIPTLOG
    fi
}

function par 
{
# Run command $1 in parallel with $3 as arguments and write output to $2, failed commands are logged to $ERRORLOG
    tolog $@
    echo $3 | tr ' ' '\n' | $CUREBIN/parallel "$1 >> $2 2>&1 || echo \"Command: \'$1\' in `pwd` failed\" >> $ERRORLOG "
}

function parX 
{
# Run command $1 in parallel with $2 as arguments and write output to $3, failed commands are logged to $ERRORLOG
    tolog $@
    echo $4 | tr ' ' '\n' | $CUREBIN/parallel -X "$1 || echo \"Command: \'$2\' in `pwd` failed\" >> $ERRORLOG " >> $3 2>&1
}

#par "$FT/meanfits --maxmem 1024 -s -m -k 2.8 -o masterbias_s.fits s{}" "`ls`" logf err

function mkbias
{
    # subtract off the overscan and averge the bias frames
    # result in $1

    # First create lists with binned and unbinned filenames
    SFILES=""
    LFILES=""
    pushd $BIAS &>/dev/null
    for f in ${NAMEBASE}*.fits; do
        isbinned $f
        if [ $? -eq 0 ]; then
            LFILES="$LFILES $f"
        else
            SFILES="$SFILES $f"
        fi
    done
    
    if [ "X$SFILES" != "X" ]; then # Binned files
        par "$FT/subtractfits -s -a -k 2.8 -o $BIASSEC_S -t -z {}" $BIASLOG "$SFILES" 
        parX "$FT/meanfits --maxmem 1024 -s -m -k 2.8 -o masterbias_s.fits s{}" meanfits $BIASLOG "$SFILES" 
        mv masterbias_s.fits ../
	mv e.masterbias_s.fits ../
    fi
    if [ "X$LFILES" != "X" ]; then # Binned files
        par "$FT/subtractfits -s -a -k 2.8 -o $BIASSEC_L -t -z {}" $BIASLOG "$LFILES" 
        parX "$FT/meanfits --maxmem 1024 -s -m -k 2.8 -o masterbias_l.fits s{}" meanfits $BIASLOG "$LFILES" 
        mv masterbias_l.fits ../
	mv e.masterbias_l.fits ../
    fi
    popd &>/dev/null
    if [ "X$DEBUG" == "Xfalse" ]; then
        rm -rf $BIAS/s${NAMEBASE}*.fits $BIAS/e.s${NAMEBASE}*.fits    
    fi
}

function subbias
{
    # subtract off the overscan and the bias frame and trim the image
    # subtract bias $1 off of files in directory $2

    # First create lists with binned and unbinned filenames
    SFILES=""
    LFILES=""
    pushd $1 &>/dev/null

    list=`ls ${NAMEBASE}*.fits 2> /dev/null`

    if [ "X$list" != "X" ]; then

        for f in ${NAMEBASE}*.fits; do
            isbinned $f
            if [ $? -eq 0 ]; then
                LFILES="$LFILES $f"
            else
                SFILES="$SFILES $f"
            fi
        done
        
        if [ "X$SFILES" != "X" ]; then # Binned files
#            echo par "$FT/subtractfits -s -a -k 2.8 -o $BIASSEC_S -t -z -f ../../../masterbias_s.fits {}" $FLATLOG "$SFILES"
            par "$FT/subtractfits -s -a -k 2.8 -o $BIASSEC_S -t -z -f ../../../masterbias_s.fits {}" $FLATLOG "$SFILES"
            par "$FT/extractfits -r $TRIMSEC_S s{}" $FLATLOG "$SFILES"
        fi
        if [ "X$LFILES" != "X" ]; then # Unbinned files
            par "$FT/subtractfits -s -a -k 2.8 -o $BIASSEC_L -t -z -f ../../../masterbias_l.fits {}" $FLATLOG "$LFILES"
            par "$FT/extractfits -r $TRIMSEC_L s{}" $FLATLOG "$LFILES"
            par "$FT/binfits -b 2,1 es{}" $FLATLOG "$LFILES"
            par "mv bes{} es{}" $FLATLOG "$LFILES"
            par "mv e.bes{} e.es{}" $FLATLOG "$LFILES"
        fi

        if [ "X$DEBUG" == "Xfalse" ]; then
            rm -f s${NAMEBASE}*.fits e.s${NAMEBASE}*.fits
        fi
    fi

    popd &>/dev/null

}

function addpnoise
{
    # add photon noise to the error frames in $1
    list=`ls $1/es${NAMEBASE}*.fits 2> /dev/null`

    if [ "X$list" != "X" ]; then
        par "$FT/addphotonnoise --gain_key GAIN1 {}" $FLATLOG "$list"
        if [ "X$DEBUG" == "Xfalse" ]; then
            rm -f $1/es${NAMEBASE}*.fits $1/e.es${NAMEBASE}*.fits
        fi
    fi
}

#function bin2x1
#{
#    # bin images 2x1, overwrite originals
#    pushd $1
#    for f in ${NAMEBASE}*.fits; do
#	size=`ls -lL $f | tr -s ' ' | cut -d ' ' -f 5`
#	if [ $size -eq 8527680 ]; then 
#	    echo Binning `ls -lL $f`
#	    $FT/binfits -b 2,1 $f
#	    mv b$f $f
#	fi
#    done
#    popd
#}

function masterarc
{
    # build master arc lamp of files in directory $2
    # result in $1
    list=`ls $2/pes${NAMEBASE}*.fits 2> /dev/null`
    
    echo "Making masterarc in $2"

    if [ "X$list" != "X" ]; then
        parX "$FT/meanfits --maxmem 1024 -s -t -m -k 2.8 -o $1 {}" masterarc $FLATLOG "$list"
        if [ "X$DEBUG" == "Xfalse" ]; then
            rm -f $2/pes${NAMEBASE}*.fits $2/e.pes${NAMEBASE}*.fits
        fi
    else
    	echo [masterarc] Warning: No files to process >> $FLATLOG
    fi
}

function mastertrace
{
    # build master trace of files in directory $2
    # result in $1
    list=`ls $2/pes${NAMEBASE}*.fits 2> /dev/null`
    
    if [ "X$list" != "X" ]; then
        parX "$FT/flatfits -r 400:550,1042:1144 -k 2.8 --norm_kappa 10 -o $1 {}" mastertrace $FLATLOG "$list"
        if [ "X$DEBUG" == "Xfalse" ]; then
            rm -f $2/pes${NAMEBASE}*.fits $2/e.pes${NAMEBASE}*.fits 
        fi
    else
    	echo [mastertrace] Warning: No files to process >> $FLATLOG
    fi
}


function deform_master
{
    # build master distortion from mastertrace $2 and masterarc $1
    # log in $3
    if [ -f $2 ]; then
        if [ -f $1 ]; then
            parX "$FT2/deformer3 $H2 $H3 -a $1 -l $LINES $DEFORM1_OPTS {}" deformer3 $3 "$2"
        else
        	echo "[deform]: ERROR: $1 not found!" >> $3
        fi
    else
        echo "[deform]: ERROR: $2 not found!" >> $3
    fi
}

function deform
{
    # build distortion correction from distortion $2 and fits in directory $3
    # trace from file $4
    # result in directory $1
    # output in $5
    # Use the fiber model file given as fallback if S/N too low to re-calculate
    list=`ls $3/pes${NAMEBASE}*.fits 2> /dev/null`
    
    if [ "X$list" != "X" ]; then
        par "$FT2/deformer3 $H2 $H3 -d $2 -o $1 -f $4 $DEFORM2_OPTS {}" $5 "$list"
    fi
}

function subtractsky
{
    # subtract skyfrom using distortion and traces in $1
    # from fits in directory $2
    list=`cd $2 && echo pes${NAMEBASE}*.fits | sed -e 's/\.fits//g'`

    if [ "$list" != "pes${NAMEBASE}*" ]; then
        par "$FT2/subtractsky -f $1/{}.fmod -d $1/{}.dist $SUBSKY_OPTS $2/{}.fits" $SKYLOG "$list"
    fi

}

function makebias
{
    mkbias 
}

function makeflat
{
    for D in $MORN $EVE; do

        subbias      $COMP/$D
        subbias      $FLAT/$D
        subbias      $SCI/$D
        addpnoise    $COMP/$D
        addpnoise    $FLAT/$D
        addpnoise    $SCI/$D
        masterarc    ${D}"_masterarc.fits" $COMP/$D
        mastertrace  ${D}"_mastertrace.fits" $FLAT/$D
        
    done
}


function makedeform
{        
    deform_master	 $MORN"_masterarc.fits" $MORN"_mastertrace.fits" $DEFORM1LOG
    deform_master    $EVE"_masterarc.fits" $EVE"_mastertrace.fits" $DEFORM2LOG
    wait
    deform           $FLAT/$MORN $MORN"_mastertrace.dist" $SCI/$MORN $MORN"_mastertrace.fmod" $DEFORM1LOG
    deform           $FLAT/$EVE $EVE"_mastertrace.dist" $SCI/$EVE $EVE"_mastertrace.fmod" $DEFORM2LOG
}

function makeskysub
{
    pushd object &>/dev/null
        
    subtractsky   ../$FLAT/$MORN $MORN 
    subtractsky   ../$FLAT/$EVE $EVE 

    popd &>/dev/null
}

if [ $MAKE_BIAS -eq 1 ]; then

    pushd $DATADIR &> /dev/null

    echo "Making bias in $DATADIR/$BIAS..."

    cat /dev/null > $BIASLOG
    makebias 

    popd &>/dev/null
fi
	
for d in $SUBDIRS; do

    pushd $d &> /dev/null
    
    if [ $MAKE_FLAT -eq 1 ]; then
        echo "Making flat in $d..."
        makeflat
    fi

    if [ $MAKE_DEFORM -eq 1 ]; then
        echo "Making deform in $d..."
        makedeform 
    fi

    if [ $MAKE_SKYSUB -eq 1 ]; then
        echo "Subtracting sky in $d..."
        makeskysub 
    fi
    popd &>/dev/null

done

