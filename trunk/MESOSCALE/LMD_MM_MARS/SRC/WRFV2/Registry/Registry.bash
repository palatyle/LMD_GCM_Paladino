#! /bin/bash
#
# This is a script to generate FORTRAN code from registry
#
# A. Spiga - 17/03/08
#


#
#
#
echo AUTOMATIC FILE GENERATION
grep '#SAVEMARS2' Registry.EM | grep 'rh' > fromregistry2
grep '#SAVEMARS3' Registry.EM | grep 'rh' > fromregistry3
yeah2=$(awk 'END {print NR}' fromregistry2)
yeah3=$(awk 'END {print NR}' fromregistry3)
echo save ${yeah2} 2D arrays
echo save ${yeah3} 3D arrays


#
#
#
\rm wrf_output_2d.h 2> /dev/null
touch wrf_output_2d.h

echo '!-------------------' >> test_include
echo '! INCLUDE 'wrf_output_2d.h'' >> test_include
echo '!-------------------' >> test_include
if [ ${yeah2} -ne 0 ]
then
awk '{print "      INTEGER, PARAMETER :: ind_"$3" = " NR}' fromregistry2 >> test_include
awk 'END {print "      INTEGER, PARAMETER :: n2d = " NR}' fromregistry2 >> test_include
mv test_include wrf_output_2d.h
else
echo '           INTEGER, PARAMETER :: n2d = 1' >> test_include
mv test_include wrf_output_2d.h
fi

\rm wrf_output_3d.h 2> /dev/null
touch wrf_output_3d.h

echo '!-------------------' >> test_include
echo '! INCLUDE 'wrf_output_3d.h'' >> test_include
echo '!-------------------' >> test_include
if [ ${yeah3} -ne 0 ]
then
awk '{print "      INTEGER, PARAMETER :: ind_"$3" = " NR}' fromregistry3 >> test_include
awk 'END {print "      INTEGER, PARAMETER :: n3d = " NR}' fromregistry3 >> test_include
mv test_include wrf_output_3d.h
else
echo '           INTEGER, PARAMETER :: n3d = 1' >> test_include
mv test_include wrf_output_3d.h
fi


#
#
#
\rm fill_save_2d.F90 2> /dev/null
#if [ ${yeah2} -ne 0 ]
#then
echo 'SUBROUTINE fill_save_2d( &' >> test_subroutine
awk '{print $NF ",&"}' fromregistry2 >> test_subroutine
echo 'ngrid,&' >> test_subroutine
echo 'output_tab2d)' >> test_subroutine
echo '   ' >> test_subroutine
echo 'IMPLICIT NONE' >> test_subroutine
echo '   ' >> test_subroutine
echo 'include "wrf_output_2d.h"' >> test_subroutine
#echo 'include "dimension.h"' >> test_subroutine
echo '   ' >> test_subroutine
echo 'INTEGER :: ngrid' >> test_subroutine
awk '{print "REAL, DIMENSION(ngrid) :: " $NF}' fromregistry2 >> test_subroutine
echo 'REAL, DIMENSION(ngrid,n2d) :: output_tab2d' >> test_subroutine
echo '   ' >> test_subroutine
awk '{print "output_tab2d(:,ind_"$3")="$NF"(:)"}' fromregistry2 >> test_subroutine
echo '   ' >> test_subroutine
echo 'END SUBROUTINE fill_save_2d' >> test_subroutine
mv test_subroutine fill_save_2d.F90
#fi

\rm fill_save_3d.F90 2> /dev/null
#if [ ${yeah3} -ne 0 ]
#then
echo 'SUBROUTINE fill_save_3d( &' >> test_subroutine
awk '{print $NF ",&"}' fromregistry3 >> test_subroutine
echo 'ngrid,&' >> test_subroutine
echo 'nlayer,&' >> test_subroutine
echo 'output_tab3d)' >> test_subroutine
echo '   ' >> test_subroutine
echo 'IMPLICIT NONE' >> test_subroutine
echo '   ' >> test_subroutine
echo 'include "wrf_output_3d.h"' >> test_subroutine
#echo 'include "dimension.h"' >> test_subroutine
echo '   ' >> test_subroutine
echo 'INTEGER :: ngrid' >> test_subroutine
echo 'INTEGER :: nlayer' >> test_subroutine
awk '{print "REAL, DIMENSION(ngrid,nlayer) :: " $NF}' fromregistry3 >> test_subroutine
echo 'REAL, DIMENSION(ngrid,nlayer,n3d) :: output_tab3d' >> test_subroutine
echo '   ' >> test_subroutine
awk '{print "output_tab3d(:,:,ind_"$3")="$NF"(:,:)"}' fromregistry3 >> test_subroutine
echo '   ' >> test_subroutine
echo 'END SUBROUTINE fill_save_3d' >> test_subroutine
mv test_subroutine fill_save_3d.F90
#fi


#
#
#
\rm fill_save.inc 2> /dev/null

#if [ ${yeah2} -ne 0 ]
#then
echo '      CALL fill_save_2d(' >> fill_save.inc
awk '{print "     . "$NF ","}' fromregistry2 >> fill_save.inc
echo '     . ngrid,' >> fill_save.inc
echo '     . output_tab2d)' >> fill_save.inc
#fi

#if [ ${yeah3} -ne 0 ]
#then
echo '      CALL fill_save_3d(' >> fill_save.inc
awk '{print "     . "$NF ","}' fromregistry3 >> fill_save.inc
echo '     . ngrid,' >> fill_save.inc
echo '     . nlayer,' >> fill_save.inc
echo '     . output_tab3d)' >> fill_save.inc
#fi



#
#
#
\rm module_lmd_driver_output1.inc 2> /dev/null
touch module_lmd_driver_output1.inc
\rm module_lmd_driver_output2.inc 2> /dev/null
touch module_lmd_driver_output2.inc
\rm module_lmd_driver_output3.inc 2> /dev/null
touch module_lmd_driver_output3.inc

echo 'REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &' >> module_lmd_driver_output2.inc

if [ ${yeah2} -ne 0 ]
then
awk '{print "       " $3 ",&"}' fromregistry2 >> module_lmd_driver_output1.inc
#more module_lmd_driver_output1.inc >> module_lmd_driver_output2.inc
cat module_lmd_driver_output1.inc >> module_lmd_driver_output2.inc
awk '{print $3"(i,j) = output_tab2d(subs,ind_"$3")"}' fromregistry2 >> module_lmd_driver_output3.inc
fi

echo '   PSFC,TSK' >> module_lmd_driver_output2.inc
echo 'REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT)  :: &' >> module_lmd_driver_output2.inc

if [ ${yeah3} -ne 0 ]
then
awk '{print "        " $3 ",&"}' fromregistry3 >> module_lmd_driver_output1.inc
awk '{print "        " $3 ",&"}' fromregistry3 >> module_lmd_driver_output2.inc
echo 'DO k = kps,kpe' >> module_lmd_driver_output3.inc
awk '{print $3"(i,k,j) = output_tab3d(subs,k,ind_"$3")"}' fromregistry3 >> module_lmd_driver_output3.inc
echo 'ENDDO' >> module_lmd_driver_output3.inc
fi

echo '   RTHBLTEN,RUBLTEN,RVBLTEN' >> module_lmd_driver_output2.inc



#
#
#
\rm module_lmd_driver_output4.inc 2> /dev/null
touch module_lmd_driver_output4.inc
awk '{print "  &    ,"$3"=grid%"$3"   &"}' fromregistry2 >> module_lmd_driver_output4.inc
awk '{print "  &    ,"$3"=grid%"$3"   &"}' fromregistry3 >> module_lmd_driver_output4.inc



#
#
#
\rm fromregistry2
\rm fromregistry3
mv module_lmd_driver_output?.inc ../inc/
\rm ../mars_lmd/libf/phymars/fill_save_?d.F90 2> /dev/null
\rm ../mars_lmd/libf/phymars/wrf_output_?d.h 2> /dev/null
mv fill_save.inc fill_save_?d.F90 wrf_output_?d.h ../mars_lmd/libf/phymars/

#
#
#
echo END of AUTOMATIC FILE GENERATION
echo I guess you modified 'Registry.EM' ... 
echo ... so remove the file 'Registry' and recompile the model
