#!/bin/bash
# script to download and install the latest version of IOIPSL on Gnome
#

#0. Preliminary stuff 

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# add a "gnome" configuration to AA_make.gdef
echo "#-Q- gnome  #- Global definitions for Idataplex (gnome) at UPMC, ifort" >> AA_make.gdef
echo "#-Q- gnome  M_K = make" >> AA_make.gdef
echo "#-Q- gnome  P_C = cpp" >> AA_make.gdef
echo '#-Q- gnome  P_O = -P -C $(P_P)' >> AA_make.gdef
echo "#-Q- gnome  F_C = ifort -mcmodel=large -shared-intel -c" >> AA_make.gdef
echo "#-Q- gnome  #-D- MD    F_D = -g" >> AA_make.gdef
echo "#-Q- gnome  #-D- MN    F_D =" >> AA_make.gdef
echo "#-Q- gnome  #-P- I4R4  F_P = -integer-size 32" >> AA_make.gdef
echo "#-Q- gnome  #-P- I4R8  F_P = -integer-size 32 -real-size 64" >> AA_make.gdef
echo "#-Q- gnome  #-P- I8R8  F_P = -integer-size 64 -real-size 64" >> AA_make.gdef
echo '#-Q- gnome  F_O = -O $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)' >> AA_make.gdef
echo "#-Q- gnome  F_L = ifort" >> AA_make.gdef
echo "#-Q- gnome  M_M = 0" >> AA_make.gdef
echo "#-Q- gnome  L_X = 0" >> AA_make.gdef
echo "#-Q- gnome  L_O =" >> AA_make.gdef
echo "#-Q- gnome  A_C = ar -r" >> AA_make.gdef
echo "#-Q- gnome  A_G = ar -x" >> AA_make.gdef
echo "#-Q- gnome  C_C = icc -c" >> AA_make.gdef
echo "#-Q- gnome  C_O =" >> AA_make.gdef
echo "#-Q- gnome  C_L = icc" >> AA_make.gdef
echo "#-Q- gnome  #-" >> AA_make.gdef
echo "#-Q- gnome  NCDF_INC = /usr/local/include" >> AA_make.gdef
echo "#-Q- gnome  NCDF_LIB = -L/usr/local/lib -lnetcdf" >> AA_make.gdef
echo "#-Q- gnome  #-" >> AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t gnome -p I4R8

## 3. build ioipsl:
cd ../modeles/IOIPSL/src
make
## Compile the rebuild tool:
cd ../tools
make

if [[ -f ${whereami}/modipsl/lib/libioipsl.a ]] 
  then
  echo "OK: ioipsl library is in ${whereami}/modipsl/lib"
else
  echo "Something went wrong..."
fi
