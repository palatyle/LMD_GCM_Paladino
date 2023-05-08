#! /bin/bash

########################
# make a minimum startbase 
# with what is 
# in start_archive.nc file for clim scenario
# downloaded from the server (reference MCD 5.2)
########################

scen="clim" # the one set in callphys.def
#scen="MY29"

namefolder="startbase_"$scen"_mini"

\rm -rf $namefolder
mkdir $namefolder

\rm start_archive_64x48x49_$scen.nc
wget -q http://www.lmd.jussieu.fr/~lmdz/planets/mars/starts/start_archive_64x48x49_$scen.nc
ln -sf start_archive_64x48x49_$scen.nc start_archive.nc

echo -e "0\n \n0\n \n" | newstart.e > $namefolder/log0 2> $namefolder/log0
mv restart.nc $namefolder/start0.nc
mv restartfi.nc $namefolder/startfi0.nc

echo -e "0\n \n193\n \n" | newstart.e > $namefolder/log193 2> $namefolder/log193
mv restart.nc $namefolder/start193.nc
mv restartfi.nc $namefolder/startfi193.nc

#echo -e "0\n \n372\n \n" | newstart.e > $namefolder/log372 2> $namefolder/log372
#mv restart.nc $namefolder/start372.nc
#mv restartfi.nc $namefolder/startfi372.nc
#
#echo -e "0\n \n515\n \n" | newstart.e > $namefolder/log515 2> $namefolder/log515
#mv restart.nc $namefolder/start515.nc
#mv restartfi.nc $namefolder/startfi515.nc

