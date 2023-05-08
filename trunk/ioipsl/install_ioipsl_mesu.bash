#!/bin/bash
# script to download and install the latest version of IOIPSL on mesu
#

#0. Preliminary stuff 
module purge
module load intel-compilers-16/16.0.2.181
module load intel-cmkl-16/16.0.2
module load mpt/2.11
module load netcdf/4.3.0

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# add a "mesu" configuration to AA_make.gdef
echo "#-Q- mesu  #- Global definitions for Idataplex (mesu) at UPMC, ifort" >> AA_make.gdef
echo "#-Q- mesu  M_K = make" >> AA_make.gdef
echo "#-Q- mesu  P_C = cpp" >> AA_make.gdef
echo '#-Q- mesu  P_O = -P -C $(P_P)' >> AA_make.gdef
echo "#-Q- mesu  F_C = ifort -mcmodel=medium -shared-intel -c" >> AA_make.gdef
echo "#-Q- mesu  #-D- MD    F_D = -g" >> AA_make.gdef
echo "#-Q- mesu  #-D- MN    F_D =" >> AA_make.gdef
echo "#-Q- mesu  #-P- I4R4  F_P = -integer-size 32" >> AA_make.gdef
echo "#-Q- mesu  #-P- I4R8  F_P = -integer-size 32 -real-size 64" >> AA_make.gdef
echo "#-Q- mesu  #-P- I8R8  F_P = -integer-size 64 -real-size 64" >> AA_make.gdef
echo '#-Q- mesu  F_O = -O $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)' >> AA_make.gdef
echo "#-Q- mesu  F_L = ifort" >> AA_make.gdef
echo "#-Q- mesu  M_M = 0" >> AA_make.gdef
echo "#-Q- mesu  L_X = 0" >> AA_make.gdef
echo "#-Q- mesu  L_O =" >> AA_make.gdef
echo "#-Q- mesu  A_C = ar -r" >> AA_make.gdef
echo "#-Q- mesu  A_G = ar -x" >> AA_make.gdef
echo "#-Q- mesu  C_C = icc -c" >> AA_make.gdef
echo "#-Q- mesu  C_O =" >> AA_make.gdef
echo "#-Q- mesu  C_L = icc" >> AA_make.gdef
echo "#-Q- mesu  #-" >> AA_make.gdef
echo "#-Q- mesu  NCDF_INC = ${NETCDF_FORT_ROOT}/include" >> AA_make.gdef
echo "#-Q- mesu  NCDF_LIB = -L${NETCDF_FORT_ROOT}/lib -lnetcdff" >> AA_make.gdef
echo "#-Q- mesu  #-" >> AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t mesu -p I4R8

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
