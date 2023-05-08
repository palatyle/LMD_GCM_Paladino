#!/bin/bash
# script to download and install the latest version of IOIPSL on occigen
#

#0. Preliminary stuff 
module purge
#module load intel/15.6.233
#module load intelmpi/5.1.3.258
#module load hdf5/1.8.14
#module load netcdf/4.3.3-rc2_fortran-4.4.1
module load intel/17.0
module load intelmpi/2017.0.098
module load hdf5/1.8.17
module load netcdf/4.4.0_fortran-4.4.2

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# add a "occigen" configuration to AA_make.gdef
echo "#-Q- occigen  #- Global definitions for Occigen at CINES, ifort" >> AA_make.gdef
echo "#-Q- occigen  M_K = make" >> AA_make.gdef
echo "#-Q- occigen  P_C = cpp" >> AA_make.gdef
echo '#-Q- occigen  P_O = -P -C $(P_P)' >> AA_make.gdef
echo "#-Q- occigen  F_C = ifort -mcmodel=medium -shared-intel -c" >> AA_make.gdef
echo "#-Q- occigen  #-D- MD    F_D = -g" >> AA_make.gdef
echo "#-Q- occigen  #-D- MN    F_D =" >> AA_make.gdef
echo "#-Q- occigen  #-P- I4R4  F_P = -integer-size 32" >> AA_make.gdef
echo "#-Q- occigen  #-P- I4R8  F_P = -integer-size 32 -real-size 64" >> AA_make.gdef
echo "#-Q- occigen  #-P- I8R8  F_P = -integer-size 64 -real-size 64" >> AA_make.gdef
echo '#-Q- occigen  F_O = -O $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)' >> AA_make.gdef
echo "#-Q- occigen  F_L = ifort" >> AA_make.gdef
echo "#-Q- occigen  M_M = 0" >> AA_make.gdef
echo "#-Q- occigen  L_X = 0" >> AA_make.gdef
echo "#-Q- occigen  L_O =" >> AA_make.gdef
echo "#-Q- occigen  A_C = ar -r" >> AA_make.gdef
echo "#-Q- occigen  A_G = ar -x" >> AA_make.gdef
echo "#-Q- occigen  C_C = icc -c" >> AA_make.gdef
echo "#-Q- occigen  C_O =" >> AA_make.gdef
echo "#-Q- occigen  C_L = icc" >> AA_make.gdef
echo "#-Q- occigen  #-" >> AA_make.gdef
echo "#-Q- occigen  NCDF_INC = ${NETCDFF_INCDIR}" >> AA_make.gdef
echo "#-Q- occigen  NCDF_LIB = -L${NETCDFF_LIBDIR} -lnetcdff" >> AA_make.gdef
echo "#-Q- occigen  #-" >> AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t occigen -p I4R8

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
