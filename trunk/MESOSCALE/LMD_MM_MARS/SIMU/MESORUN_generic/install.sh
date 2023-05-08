#! /bin/bash

## INSTALLER: MESOSCALE + PHYSTD as LIB
## A. Spiga -- 03/2020

############
## SETTINGS
############
kind=lesmpifort_64

#############
## LOCATIONS
#############
here=$PWD
mm=$here/code/MESOSCALE/LMD_MM_MARS
fd=$mm/generic_lmd_new_les_mpifort_64/WRFV2
src=$mm/SRC/
mod=$src/DEV/simpler_compile_LES_phys/

#############################
## RECORD WHAT IS BEING DONE
#############################
exec 1> install.log
exec 2> install.loge

#################
## DOWNLOAD CODE
#################
## for physics, both GENERIC/MARS *and* COMMON are needed
\rm -rf $here/code
svn co -N http://svn.lmd.jussieu.fr/Planeto/trunk code
cd $here/code
svn update LMDZ.COMMON
svn update LMDZ.GENERIC
svn update LMDZ.MARS
svn update MESOSCALE
svn info

###################################
## COMPILE PHYSICS AS SEPARATE LIB
###################################
## this is without IOIPSL
## to compile with IOIPSL (e.g. for testing)
## ... use --revision=2576
## ... in install_ioipsl_ciclad.bash: svn co --revision=2576 http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
## ... ./makelmdz_fcm -t 1 -p std -b 38x36 -full -s 1 -d 25 -arch CICLADifort -cpp MESOSCALE gcm -j 8 -libphy
## ... then use specific configure.wrf
cd $here/code/LMDZ.COMMON
./makelmdz_fcm -t 1 -p std -b 38x36 -full -s 1 -d 25 -arch CICLADifort -cpp MESOSCALE gcm -io noioipsl -j 8 -libphy

############################################
## DOWNLOAD WRF AND CREATE MODIFIED VERSION
############################################
cd $src/LES/
./LMD_LES_MARS_install

##########################################
## CREATE THE SPECIFIC FOLDER FOR COMPILE
##########################################
## only use of makemeso to install the folder
## -- configure.wrf is created but this is overridden below
## -- makemeso is not compiling (this script is sort of makemeso 2.0)
cd $mm
./makemeso -c les -p generic_lmd_new -d < $here/${kind}.compile

###########################################################################
## COPY WHAT IS NECESSARY TO HAVE CONSISTENCY BETWEEN physics AND dynamics
###########################################################################
cp $mod/external/io_grib_share/wrf_io_flags.h $fd/external/io_grib_share/wrf_io_flags.h
cp $mod/external/RSL_LITE/module_dm.F $fd/external/RSL_LITE/module_dm.F

###############################################################
## PATCH THE configure.wrf FOR INTERFACE AND COMPILE THE MODEL
###############################################################
# TODO: link to env variable in the patch instead of absolute links
cd $fd
./configure < $here/${kind}.configure
patch -b < $here/${kind}.patch
./compile em_les > log_compile 2> log_error
cd $here

########################################
## GET THE STATIC DATA NECESSARY TO RUN
########################################
# NB: the latest / is important
wget -r --no-parent -nH --cut-dirs=3  https://www.lmd.jussieu.fr/~lmdz/planets/LMDZ.GENERIC/datagcm/
rm -rf robots.txt
### one file is missing from the LMDZ repo, get from HITRAN directly
cd datagcm/continuum_data/
wget https://www.hitran.org/data/CIA/N2-N2_2011.cia
cd ../../
