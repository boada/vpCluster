#!/bin/bash

BIAS="zero"
FLAT="flat"
COMP="comp"
SCI="object"
MORN="morn"
EVE="eve"
FOCUS="focus"
TEST="test"
DARK="dark"


# Check that the necessary environment variables are set
if [ x$CUREBIN == "x" ]; then
    echo "Please set CUREBIN environment variable. e.g."
    echo "bash : export CUREBIN=~/cure/bin
    echo "csh  : setenv CUREBIN "~/cure/bin"
        exit 1
fi

# Check that subtractfits is in place (we will assume that the rest
#  of the fitstools is as well if subtractfits is.
if [ ! -x $CUREBIN/headfits ]; then
    echo "Unable to find $CUREBIN/headfits, please check that CUREBIN is"
    echo "  set correctly."
    exit 1
fi



getkeyword()
{
	if [ -f $1 ]; then
    	# extract value of keyword $2 from file $1
    	#$FT/headfits $1 | grep $2 | cut -b 12-19 | sed 's/ //g'
    	$CUREBIN/headfits -r -k $2 $1
    fi
}

sort_evening_morning()
{
    # sort files into 2 directories, evening and morning, according to
    # UT header keyword
    pushd redux/$1/$2
    mkdir -p $EVE
    mkdir -p $MORN
    
	for i in *.fits; do
		tester=`getkeyword $i "UT"`
		if [[ "$tester" > "18:00:00.00" ||  "$tester" < "06:00:00.00" ]]; then
        	    ln -f -s ../../../../../raw/$1/$i $EVE/; rm $i
		fi
		if [[ "$tester" > "06:00:00.00" && "$tester"  < "18:00:00.00" ]]; then
        	    ln -f -s ../../../../../raw/$1/$i $MORN/; rm $i
		fi
  	done

    if [ $3 -eq 1 ]; then # Create links for flats and comp only

    # Checks added by JJA on 6-7-08 to check for absent eve/morn data and 
    # if absent simply copy over morn/eve data in its place
        morn_num=`ls $MORN/ | wc -l`
        eve_num=`ls $EVE/ | wc -l`
        if [[ "$morn_num" -gt "0" && "$eve_num" -lt "1" ]]; then
            for j in $MORN/* ; do ln -f -s ../$j $EVE/`basename $j`; done
        else
            if [[ "$morn_num" -lt "1" && "$eve_num" -gt "0" ]]; then
                for j in $EVE/* ; do ln -f -s ../$j $MORN/`basename $j`; done
            fi
        fi
    fi
    popd
}

if [ ! -d raw ]; then
    echo "No raw data directory found!"
    exit 1
fi

#for obsrun in `cd raw; ls -d *`; do
for obsrun in $@; do


    mkdir -p redux/${obsrun}/${BIAS}

    for day in `cd raw/${obsrun}; ls -d 20??????`; do # Loop over the month

        mkdir -p redux/${obsrun}/${day}/$FLAT
        mkdir -p redux/${obsrun}/${day}/$COMP
        mkdir -p redux/${obsrun}/${day}/$SCI
    
        for i in raw/${obsrun}/${day}/*.fits
        do
            case `getkeyword $i "IMAGETYP"` in
	        $FLAT)  ln -f -s ../../../../$i redux/${obsrun}/${day}/$FLAT/ ;;
	        $COMP)  ln -f -s ../../../../$i redux/${obsrun}/${day}/$COMP/ ;;
	        $SCI)   ln -f -s ../../../../$i redux/${obsrun}/${day}/$SCI/ ;;
           	$BIAS)  ln -f -s ../../../$i redux/${obsrun}/$BIAS ;;
           	$FOCUS) echo "Skipping focus image: " $i;;
           	$TEST)  echo "Skipping test image: " $i;;
           	$DARK)  echo "Skipping dark image: " $i;;
	        *) echo "Unknown image type" $i ;;
            esac
        done
    
# remove dome flats
        pushd redux/${obsrun}/${day}/$FLAT
        files=*.fits                                                                                                                                                                   
        if [[ $files != "*.fits" ]]; then                                                                                                                                                    
        	$CUREBIN/headfits -f -r -k OBJECT *.fits | grep dome | cut -f 1 -d " " | xargs rm
        fi
        popd
    
# remove diffuser flats
        pushd redux/${obsrun}/${day}/$SCI
        files=*.fits                                                                                                                                                                   
        if [[ $files != "*.fits" ]]; then                                                                                                                                                    
	        $CUREBIN/headfits -f -r -k OBJECT *.fits | grep flat | cut -f 1 -d " " | xargs rm
	    fi
        popd

# sort the flats and comps into evening and morning
        sort_evening_morning ${obsrun}/${day} $COMP 1
        sort_evening_morning ${obsrun}/${day} $FLAT 1 
        sort_evening_morning ${obsrun}/${day} $SCI 0
        
    done

done
