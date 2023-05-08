#!/bin/bash
# script to download and install the latest version of IOIPSL on occigen
#

#0. Preliminary stuff 
source ../arch/arch-X64_JEANZAY.env
module load subversion


whereami=`pwd -P`

# 1. Get IOIPSL (via modipsl)
svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl
cd modipsl/util

./model IOIPSL_PLUS

# 2. Set correct settings:
# set default working precision for IOIPSL:
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
