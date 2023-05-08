#!/bin/bash -x
# script to download and install the latest version of IOIPSL
# using gfortran
# You'll probably have to change paths to NetCDF library 'lib' and 'include'
# below to adapt this script to your computer.

#setfolder="/cm/shared/apps/netcdf/gcc/64/4.1.1"
#setfolder="/donnees/emlmd/netcdf64-4.0.1_gfortran/"
setfolder="/usr/local/netcdf-4.0.1/"

#0. Preliminary stuff 
# netcdf include and lib dirs:
netcdf_include=$setfolder"/include"
netcdf_lib=$setfolder"/lib"

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
export MAKE=make
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# Add a Darwin case in 'AA_make.gdef'

echo '#-Q- Darwin  #- Global definitions for gfortran on a Mac' >> AA_make.gdef
echo '#-Q- Darwin  M_K = make' >> AA_make.gdef
echo '#-Q- Darwin  P_C = cpp' >> AA_make.gdef
echo '#-Q- Darwin  P_O = -P -C -traditional $(P_P)' >> AA_make.gdef
echo '#-Q- Darwin  F_C = gfortran -c -cpp' >> AA_make.gdef
echo '#-Q- Darwin  #-D- MD    F_D = -g -Wall -fbounds-check -pedantic -finit-real=nan' >> AA_make.gdef
echo '#-Q- Darwin  #-D- MN    F_D =' >> AA_make.gdef
echo '#-Q- Darwin  #-P- I4R4  F_P =' >> AA_make.gdef
echo '#-Q- Darwin  #-P- I4R8  F_P = -fdefault-real-8' >> AA_make.gdef
echo '#-Q- Darwin  #-P- I8R8  F_P = -fdefault-integer-8 -fdefault-real-8' >> AA_make.gdef
echo '#-Q- Darwin  w_w = -O3 -funroll-all-loops $(F_D) $(F_P) -I$(MODDIR)' >> AA_make.gdef
echo '#-Q- Darwin  F_O = $(w_w) -J$(MODDIR)' >> AA_make.gdef
echo '#-Q- Darwin  F_L = gfortran' >> AA_make.gdef
echo '#-Q- Darwin  M_M = 0' >> AA_make.gdef
echo '#-Q- Darwin  L_X = 0' >> AA_make.gdef
echo '#-Q- Darwin  L_O =' >> AA_make.gdef
echo '#-Q- Darwin  A_C = ar -rs' >> AA_make.gdef
echo '#-Q- Darwin  A_G = ar -x' >> AA_make.gdef
echo '#-Q- Darwin  C_C = cc -c' >> AA_make.gdef
echo '#-Q- Darwin  C_O =' >> AA_make.gdef
echo '#-Q- Darwin  C_L = cc' >> AA_make.gdef
echo '#-Q- Darwin  #-' >> AA_make.gdef
echo "#-Q- Darwin  NCDF_INC = ${netcdf_include}" >> AA_make.gdef
echo "#-Q- Darwin  NCDF_LIB = -L${netcdf_lib} -lnetcdf -lnetcdff" >> AA_make.gdef
echo '#-Q- Darwin  #-' >> AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t Darwin -p I4R8

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

