#! /bin/bash


logfile="diff_meso_gcm.log"
wheregcm=$1

\rm $logfile
touch $logfile
#echo '***** meso_inifis.F' >> $logfile
diff ./meso_inifis.F      $1/inifis.F          >> $logfile
#echo '***** meso_callkeys.h' >> $logfile 
diff ./meso_callkeys.h    $1/callkeys.h        >> $logfile
#echo '***** meso_dustlift.F' >> $logfile
#diff ./meso_dustlift.F    $1/dustlift.F   >> $logfile
#echo '***** meso_inifis.F' >> $logfile
diff ./meso_inifis.F      $1/inifis.F          >> $logfile
#echo '***** meso_newcondens.F' >> $logfile
diff ./meso_newcondens.F  $1/newcondens.F      >> $logfile
#echo '***** meso_physiq.F' >> $logfile
diff ./meso_physiq.F      $1/physiq.F          >> $logfile
diff ./meso_testphys1d.F  $1/testphys1d.F      >> $logfile
