#!/bin/bash
# script to download and install the latest version of IOIPSL
# using pgf90
# You'll probably have to change paths to NetCDF library 'lib' and 'include'
# below to adapt this script to your computer.

#0. Preliminary stuff 
# netcdf include and lib dirs:
netcdf_include="/donnees/emlmd/netcdf64-4.0.1_pgi/include"
netcdf_lib="/donnees/emlmd/netcdf64-4.0.1_pgi/lib"

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# modify path to netcdf in 'AA_make.gdef'
cp AA_make.gdef AA_make.gdef.old
sed -e s:"linux    NCDF_INC = /distrib/local/netcdf/pgf/include/":"linux    NCDF_INC = ${netcdf_include}":1 \
    -e s:"linux    NCDF_LIB = -L/distrib/local/netcdf/pgf/lib/":"linux    NCDF_LIB = -L${netcdf_lib}":1 \
    AA_make.gdef.old > AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t linux -p I4R8

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

