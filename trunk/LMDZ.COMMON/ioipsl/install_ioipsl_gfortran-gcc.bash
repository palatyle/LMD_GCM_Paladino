#!/bin/bash
# script to download and install the latest version of IOIPSL
# using gfortran
# You'll probably have to change paths to NetCDF library 'lib' and 'include'
# below to adapt this script to your computer.

#setfolder="/cm/shared/apps/netcdf/gcc/64/4.1.1"
#setfolder="/donnees/emlmd/netcdf64-4.0.1_gfortran/"

#0. Preliminary stuff 
# netcdf include and lib dirs:
#netcdf_include=$setfolder"/include"
#netcdf_lib=$setfolder"/lib"
netcdf_include="$NETCDFINCLUDE"
netcdf_lib="$NETCDFDIR"
echo $netcdf_include
echo $netcdf_lib

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

# make all refs to ksh become refs to bash
for i in ins_m_prec model script_diff_model script_log_analyse script_recup_model
do
  sed -i -e s:'#!/bin/ksh':'#!/bin/bash':1 $i
done
./model IOIPSL_PLUS

# 2. Set correct settings:
# modify path to netcdf in 'AA_make.gdef'
cp AA_make.gdef AA_make.gdef.old
sed -e s:"gfortran  NCDF_INC = /usr/local/include":"gfortran  NCDF_INC = ${netcdf_include}":1 \
    -e s:"gfortran  NCDF_LIB = -L/usr/local/lib -lnetcdf":"gfortran  NCDF_LIB = -L${netcdf_lib} -lnetcdf -lnetcdff":1 \
    AA_make.gdef.old > AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t gfortran -p I4R8

## 3. build ioipsl:
cd ../modeles/IOIPSL/src
# make all refs to ksh become refs to bash
for i in AA_make.ldef Makefile
do
  sed -i -e s:'/bin/ksh':'/bin/bash':1 $i
done
make
## Compile the rebuild tool:
cd ../tools
for i in AA_make.ldef Makefile rebuild
do
  sed -i -e s:'/bin/ksh':'/bin/bash':1 $i
done
make

if [[ -f ${whereami}/modipsl/lib/libioipsl.a ]] 
  then
  echo "OK: ioipsl library is in ${whereami}/modipsl/lib"
else
  echo "Something went wrong..."
fi

