#!/bin/bash
# script to download and install the latest version of IOIPSL on Ciclad
#

#0. Preliminary stuff
# source the environment from the GCM arch files
source ../arch/arch-CICLADifort.env
source ../arch/arch-CICLADifort.path

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# add a "ciclad" configuration to AA_make.gdef
echo "#-Q- cicladi  #- Global definitions for ciclad at UPMC, ifort" >> AA_make.gdef
echo "#-Q- cicladi  M_K = make" >> AA_make.gdef
echo "#-Q- cicladi  P_C = cpp" >> AA_make.gdef
echo '#-Q- cicladi  P_O = -P -C $(P_P)' >> AA_make.gdef
echo "#-Q- cicladi  F_C = ifort -mcmodel=large -shared-intel -c" >> AA_make.gdef
echo "#-Q- cicladi  #-D- MD    F_D = -g" >> AA_make.gdef
echo "#-Q- cicladi  #-D- MN    F_D =" >> AA_make.gdef
echo "#-Q- cicladi  #-P- I4R4  F_P = -integer-size 32" >> AA_make.gdef
echo "#-Q- cicladi  #-P- I4R8  F_P = -integer-size 32 -real-size 64" >> AA_make.gdef
echo "#-Q- cicladi  #-P- I8R8  F_P = -integer-size 64 -real-size 64" >> AA_make.gdef
echo '#-Q- cicladi  F_O = -O $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)' >> AA_make.gdef
echo "#-Q- cicladi  F_L = ifort" >> AA_make.gdef
echo "#-Q- cicladi  M_M = 0" >> AA_make.gdef
echo "#-Q- cicladi  L_X = 0" >> AA_make.gdef
echo "#-Q- cicladi  L_O =" >> AA_make.gdef
echo "#-Q- cicladi  A_C = ar -r" >> AA_make.gdef
echo "#-Q- cicladi  A_G = ar -x" >> AA_make.gdef
echo "#-Q- cicladi  C_C = icc -c" >> AA_make.gdef
echo "#-Q- cicladi  C_O =" >> AA_make.gdef
echo "#-Q- cicladi  C_L = icc" >> AA_make.gdef
echo "#-Q- cicladi  #-" >> AA_make.gdef
echo "#-Q- cicladi  NCDF_INC = ${NETCDF_INCDIR:2}" >> AA_make.gdef
echo "#-Q- cicladi  NCDF_LIB = ${NETCDF_LIBDIR} ${NETCDF_LIB}" >> AA_make.gdef
echo "#-Q- cicladi  #-" >> AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t cicladi -p I4R8

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
