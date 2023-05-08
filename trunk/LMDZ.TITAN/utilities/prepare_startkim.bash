#!/bin/bash

################################################################################
#
# Author : Jan Vatant d'Ollone (2018)
# ------
#
# Purpose : * Modify Titan startfi.nc files to fill them with upper chemistry
# -------   fields from old comp_XX files that were used in ye old days with
#           a 32x48x55 GCM grid which implied 70 upper levels.
#
#           * You can either just fill the fields or create them if needed
#
#           * Allows you to do the same for pressure grid "preskim" field :
#            just give a vertical grid as second input file !i
#
#           * You need to run it in a folder containing comp_01->49 files !
#
#           * Requirements : NCO tool ncap2
#
# Inputs  : $1 : startfi.nc file
# ------
#          $2 (optional) : upper mid-layers pressure grid ASCII file
#          NB : It should be 1D with decreasing pressure
#
# Execution exemple : ./prepare_startkim.bash startfi.nc preskim.dat 
# -----------------
#
# Approx. run time ~ 1h30-2h
# ----------------
#
###############################################################################

ncfile="$1" # startfi.nc file

tmpfile="tmp" # temporary file


# User choice :
# * 1 Create upper chemistry dimension and fields
# * 2 Just fill them, assuming you already have them in file

while true; do
    read -p "Choose : - 1 - Create upper chemistry dimension and fields - 2 - Just fill fields assuming they pre-exist : " yn
    case $yn in
	[1]* ) echo "Ok I will create dimension upper_chemistry_layers=70 and upper chemistry fields..." ; create=true ; break ;;
	[2]* ) echo "Ok I assume you have upper chemistry fields in startfi.nc file, I will just fill them..."; create=false; break ;;
	* ) echo "Please enter 1 or 2 !";;
    esac
done

# If needed, we create the 70 levels dimension in the startfi.nc file

if [ "$create" = true ] ; then
    echo "I add to startfi a upper_chemistry_layers=70 dimension ..."
    ncap2 -s 'defdim("upper_chemistry_layers",70)' -O $ncfile
fi

# If a second argument is given it will create and fill preskim

if [ -z "$2" ] ; then
    echo " No second argument found, skipping pressure grid stuff ... "
else
    presfile="$2"
    if [ "$create" = true ] ; then
	echo "I create an empty preskim field..."
	ncap2 -s 'preskim[$upper_chemistry_layers]=0.0;preskim@title="Upper chemistry mid-layer pressure"' $ncfile
    fi

    echo "I fill preskim field with input profile ..."
    i=0
    while read -r pres
    do
	ncap2 -s "preskim($i)=$pres" -h -A $ncfile
	i=$((i+1))
    done < "$presfile"
fi

# Loop on all 44 Titan chemistry species
echo "Starting loop on all 44 Titan chemistry species ..."
for ((iq=1;iq<=44;iq++)) ; do
    
    # First calculates lines where specie name and values are written in old comp file
    qline=$(( 2 + 71*$((iq-1)) )) # Specie name
    
    qdeb=$((qline+1))              # First line of composition for this specie
    qend=$((qline+70))             # Last   "   "     "         "   "    "
    
    specie=$(sed -n "${qline}p" comp_01) # Read specie name
    specie=$(echo $specie) # Trim it
    
    echo $specie

    # Create field if needed
    if [ "$create" = true ] ; then
	echo "I create an empty upper chemistry field for", $specie
	ncap2 -s "${specie}_up[$upper_chemistry_layers,$physical_points]=0.0;${specie}_up@title='${specie} in upper atmosphere'" $ncfile
    fi
    
    # Loop on latitudes
    
    for ((ilat=1;ilat<=49;ilat++)) ; do
	
	# Deal with comp_0X filenames
	if [ $ilat -lt 10 ] ; then  
	    lat=0$ilat
	else
	    lat=$ilat
	fi

	# Extract the current specie section in latitude file
	sed -n "${qdeb},${qend}p" comp_$lat > tmp
	
	i=0
	
	echo "Filling : " $specie "-" comp_$lat

	# Deal with special mono-gridpoint for North and South Poles
	if [ $ilat -eq 1 ] ; then
	    ngrid0=0
	    ngrid1=0
	elif [ $lat -eq 49 ] ; then
	    ngrid0=1505
	    ngrid1=1505
	else
	    ngrid0=$((1+$(($ilat-2))*32))
	    ngrid1=$(($ngrid0+31))
	fi
	
	# Read tmp file and fill the grid points 
	while read -r dum1 dum2 ykim
	do
	    ncap2 -s "${specie}_up($i,$ngrid0:$ngrid1)=$ykim" -h -A $ncfile
	    i=$((i+1))
	done < "$tmpfile"

    done # Latitude loop

    rm tmp

done # Specie loop

echo "Everything's finished !"
