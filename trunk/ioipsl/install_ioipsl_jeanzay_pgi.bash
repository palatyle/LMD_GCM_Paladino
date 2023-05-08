#!/bin/bash
# script to download and install the latest version of IOIPSL on occigen
#

#0. Preliminary stuff 
source ../arch/arch-X64_JEANZAY-pgi.env
module load subversion


whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# set default working precision for IOIPSL:
sed -i -e s:"jeanzay  F_C = mpiifort -c -cpp":"jeanzay  F_C = mpif90 -c -Mpreprocess":1 \
-e s/'jeanzay  F_O = -DCPP_PARA -O3 $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR) -fp-model precise'/'jeanzay  F_O = -DCPP_PARA -fast -O3 -Munroll=c:4 $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)'/1 \
-e s:'jeanzay  NCDF_INC = ./':"jeanzay  NCDF_INC = $(nf-config --includedir)":1 \
-e s:'jeanzay  NCDF_LIB = -lnetcdff':"jeanzay  NCDF_LIB = $(nf-config --flibs)":1 \
-e s:"jeanzay  F_L = mpiifort":"jeanzay  F_L = mpif90":1 AA_make.gdef

./ins_make -t jeanzay -p I4R8

## 3. build ioipsl:
cd ../modeles/IOIPSL/src
gmake
## Compile the rebuild tool:
cd ../tools
gmake

if [[ -f ${whereami}/modipsl/lib/libioipsl.a ]] 
  then
  echo "OK: ioipsl library is in ${whereami}/modipsl/lib"
else
  echo "Something went wrong..."
fi
