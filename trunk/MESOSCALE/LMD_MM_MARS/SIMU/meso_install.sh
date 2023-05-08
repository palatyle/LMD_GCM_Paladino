#! /bin/bash

################################
# A. Spiga 09/06/2015
# Install the LMD mesoscale model
################################
# prerequisite: NETCDF
# -- NETCDF env variable
################################
## DEFAULT
## name of the folder
name="MESORUN"
## SVN version
version="HEAD"
version="1520"
################################

meso="on"
gcm=""
gcmcompile="on"
les="off"
fortcom="ifort"
while getopts "gn:hv:f:ol" options; do
  case $options in
   n ) name="${OPTARG}";;
   g ) meso="off";; 
   v ) version="${OPTARG}";;
   o ) gcm="old" ;;
   l ) gcmcompile="off";les="on" ;;
   f ) fortcom="${OPTARG}" ;;
   h ) echo "
# OPTIONS
# -n [name]        --> name of the folder to be created
# -g               --> only compile GCM (for tests)
# -v [# (or) HEAD] --> svn version 
# -o               --> old GCM+meso version
# -l               --> LES version (only new physics)
# -f               --> fortran compiler (experimental)
" ; exit ;;
  esac
done

################################
## machine on which you will compile
machine="CICLAD"
## server for sources
where_is_svn="http://svn.lmd.jussieu.fr/Planeto/trunk/"
## grid definition for GCM
dimgcm="64x48x29"
## location of static data
webrepo="http://www.lmd.jussieu.fr/~aslmd/mesoscale_model/data_static/"
## TBD: datadir:  http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir/
################################

## -----------------------------
## import settings and structure
## -----------------------------
echo "*** get structure"
rm -rf $name > /dev/null 2> /dev/null
svn -q co $where_is_svn"/MESOSCALE/LMD_MM_MARS/SIMU/MESORUN"$gcm $name
refrepo=$PWD/$name
## fill here user input to obtain independent script
case ${fortcom}$gcm in
  "ifort")    echo 1 > $refrepo/wpsin ; echo 5 > $refrepo/mesoin ; echo 1 >> $refrepo/mesoin ;;
  "ifortold") echo 1 > $refrepo/wpsin ; echo 5 > $refrepo/mesoin ; echo 4 >> $refrepo/mesoin ;
              echo 61 >> $refrepo/mesoin ; echo 61 >> $refrepo/mesoin ; echo 61 >> $refrepo/mesoin ; 
              echo 1 >> $refrepo/mesoin ; echo 1 >> $refrepo/mesoin ;;
  "gnuold")   echo 10 > $refrepo/wpsin ; echo 8 > $refrepo/mesoin ;
              echo 61 >> $refrepo/mesoin ; echo 61 >> $refrepo/mesoin ; echo 61 >> $refrepo/mesoin ;
              echo 1 >> $refrepo/mesoin ; echo 1 >> $refrepo/mesoin ;;
  *) echo "compiler not supported" ; exit ;;
esac

## ----------------
## create code repo
## ----------------
echo "*** get SVN repository"
\rm $refrepo/code
svn -q co -N $where_is_svn $refrepo/code


###################################
################################### GCM
###################################
if [[ "${gcm}" == "old" ]]
then

 if [[ "${gcmcompile}" == "on" ]]
 then
   log=$refrepo/code/MESOSCALE/LMDZ.MARS/logcompile_gcm
   echo "*** get and compile GCM code version "$version
   cd $refrepo/code
   svn update -r $version MESOSCALE > /dev/null
   cd $refrepo/code/MESOSCALE/LMDZ.MARS
   ln -sf makegcm_$fortcom makegcm
   ./compile > $log 2> $log
 fi
 
else

 ## ------------
 ## get GCM code
 ## ------------
 echo "*** get GCM code version "$version
 cd $refrepo/code
 # at least get Mars physics (always needed)
 svn -q update -r $version LMDZ.MARS 

## START compiling GCM PART
if [[ "${gcmcompile}" == "on" ]]
then

 svn -q update -r $version LMDZ.COMMON
 cd $refrepo/code/LMDZ.COMMON
 svn -q co http://forge.ipsl.jussieu.fr/fcm/svn/PATCHED/FCM_V1.2
 ln -sf FCM_V1.2/bin/fcm .

 ## --------------
 ## compile IOIPSL
 ## --------------
 log=$refrepo/code/logcompile_ioipsl
 echo "*** compile IOIPSL: check progress in "$log
 rm -rf $log ; touch $log
 cd $refrepo/code/LMDZ.COMMON/ioipsl
 if [[ "${machine}" == "CICLAD" ]]
 then
   ./install_ioipsl_ciclad.bash > $log 2> $log
 else
   ./install_ioipsl_$fortcom".bash" > $log 2> $log
 fi

 ## -----------
 ## compile GCM
 ## -----------
 log=$refrepo/code/logcompile_gcm
 echo "*** compile GCM: check progress in "$log
 rm -rf $log ; touch $log
 # make a re-usable command
 echo "#! /bin/bash" > $refrepo/compile_gcm.sh
 echo "cd $refrepo/code/LMDZ.COMMON" >> $refrepo/compile_gcm.sh
 echo "./makelmdz_fcm -cpp MESOINI -j 8 -s 2 -d $dimgcm -arch $machine$fortcom -parallel mpi -p mars gcm" >> $refrepo/compile_gcm.sh
 echo "./makelmdz_fcm              -j 8 -s 2 -d $dimgcm -arch $machine$fortcom               -p mars newstart" >> $refrepo/compile_gcm.sh
 echo "cd $refrepo/gcm ; \rm gcm.e ; ln -sf $refrepo/code/LMDZ.COMMON/bin/gcm_${dimgcm}_phymars_para.e gcm.e" >> $refrepo/compile_gcm.sh
 echo "cd $refrepo/gcm/newstart ; \rm newstart.e ; ln -sf $refrepo/code/LMDZ.COMMON/bin/newstart_${dimgcm}_phymars_seq.e newstart.e" >> $refrepo/compile_gcm.sh 
 chmod 755 $refrepo/compile_gcm.sh
 # now execute command
 $refrepo/compile_gcm.sh > $log 2> $log

 ## ------------------------
 ## make a minimal startbase
 ## ------------------------
 echo "*** make a minimal startbase"
 cd $refrepo/gcm/newstart
 ./mini_startbase.sh

fi

fi

###################################
################################### MESO
###################################

## START MESOSCALE PART
if [[ "${meso}" == "on" ]]
then

  ###
  if [[ "${gcm}" == "old" ]]
  then
    option=""
  else
    option="-p mars_lmd_new"
  fi
  ###

## ----------------------
## get and make mesoscale
## ----------------------
echo "*** get and compile mesoscale version "$version
cd $refrepo/code
svn update -r $version MESOSCALE > /dev/null
#
if [[ "${les}" == "on" ]]
then
  echo "*** LES LES LES LES"
  cd $refrepo/code/MESOSCALE/LMD_MM_MARS/SRC/LES/
  ./LMD_LES_MARS_install > /dev/null
  option=$option" -c les"
fi
#
cd $refrepo/code/MESOSCALE/LMD_MM_MARS
ls $refrepo/mesoin
if [[ "$?" == 0 ]] ; then
  ./makemeso $option < $refrepo/mesoin
else
  ./makemeso $option
fi
rm -rf $refrepo/mesoin

if [[ "${les}" == "off" ]]
then
## -------------------------------
## make ini&bdy tools in mesoscale
## -------------------------------
echo "*** compile initialization tools"
cd $refrepo/code_compiled
ln -sf $refrepo/code/MESOSCALE/LMD_MM_MARS/SRC/SCRIPTS/prepare_ini .
./prepare_ini > /dev/null
##
cd $refrepo/code_compiled/PREP_MARS
./compile"_"$fortcom
##
cd $refrepo/code_compiled/WPS
ls $refrepo/wpsin
if [[ "$?" == 0 ]] ; then
  ./configure < $refrepo/wpsin > /dev/null 2> /dev/null
else
  ./configure
fi
rm -rf $refrepo/wpsin
rm -rf logcompile
./compile > logcompile 2>&1

## ------------------
## import static data
## ------------------
echo "*** get static data"
rm -rf $refrepo/data_static
svn co -q $where_is_svn/MESOSCALE/LMD_MM_MARS/WPS_GEOG $refrepo/data_static
cd $refrepo/data_static
rm -rf logdown
wget $webrepo"/albedo_TES.tar.gz" -a logdown
wget $webrepo"/mola_topo64.tar.gz" -a logdown
wget $webrepo"/thermal_TES.tar.gz" -a logdown
for fff in *.tar.gz; do
  tar xzvf $fff > /dev/null
  rm -rf $fff
done
fi

## ------------------------
## get and compile postproc
## ------------------------
echo "*** get and compile post-processing tool"
cd $refrepo
svn co -q https://github.com/aymeric-spiga/api/trunk postproc
cd $refrepo/postproc
./compile

fi
## END MESOSCALE PART

## -----
## check
## -----
echo "*** CHECKLIST:"
ls -lL $refrepo/gcm/gcm.e
ls -lL $refrepo/geogrid/geogrid.exe
ls -lL $refrepo/metgrid/metgrid.exe
ls -lL $refrepo/prep/readmeteo.exe
ls -lL $refrepo/data_static/albedo_TES
ls -lL $refrepo/real.exe
ls -lL $refrepo/wrf.exe
ls -lL $refrepo/postproc/api

