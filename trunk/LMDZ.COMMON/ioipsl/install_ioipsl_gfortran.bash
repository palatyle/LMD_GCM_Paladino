#!/bin/bash
# script to download and install the latest version of IOIPSL
# using gfortran
# You'll probably have to change paths to NetCDF library 'lib' and 'include'
# below to adapt this script to your computer.

#setfolder="/planeto/emlmd/netcdf64-4.0.1_gfortran4.4.7"
#setfolder="/donnees/emlmd/netcdf64-4.0.1_gfortran/"
setfolder="/apps/local/easybuild/software/netCDF-Fortran/4.4.2-gmvolf-5.5.4"

#0. Preliminary stuff 
# netcdf include and lib dirs:
netcdf_include=$setfolder"/include"
netcdf_lib=$setfolder"/lib"

whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

# make all refs to ksh become refs to bash
for i in ins_m_prec model script_diff_model script_log_analyse script_recup_model
do
  sed -e s:'#!/bin/ksh':'#!/bin/bash':1 $i > tmp
  mv -f tmp $i
done
chmod u=rwx model
chmod u=rwx ins_m_prec
./model IOIPSL_PLUS

# 2. Set correct settings:
# modify path to netcdf in 'AA_make.gdef'
cp AA_make.gdef AA_make.gdef.old
sed -e s:"gfortran  NCDF_INC = /usr/local/include":"gfortran  NCDF_INC = ${netcdf_include}":1 \
    -e s:"gfortran  NCDF_LIB = -L/usr/local/lib":"gfortran  NCDF_LIB = -L${netcdf_lib}":1 \
    AA_make.gdef.old > AA_make.gdef

# set default working precision for IOIPSL:
./ins_make -t gfortran -p I4R8

## 3. build ioipsl:
cd ../modeles/IOIPSL/src
# make all refs to ksh become refs to bash
for i in AA_make.ldef Makefile
do
  sed -e s:'/bin/ksh':'/bin/bash':1 $i > tmp
  mv -f tmp $i
done
make
## Compile the rebuild tool:
cd ../tools
# make all refs to ksh become refs to bash
for i in AA_make.ldef Makefile rebuild
do
  sed -e s:'/bin/ksh':'/bin/bash':1 $i > tmp
  mv -f tmp $i
done
chmod u=rwx rebuild
make

if [[ -f ${whereami}/modipsl/lib/libioipsl.a ]] 
  then
  echo "OK: ioipsl library is in ${whereami}/modipsl/lib"
else
  echo "Something went wrong..."
fi

