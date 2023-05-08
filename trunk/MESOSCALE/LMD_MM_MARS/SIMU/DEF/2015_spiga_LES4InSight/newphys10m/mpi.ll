# Nom du travail LoadLeveler
# @ job_name= LES
# Fichier de sortie standard du travail
# @ output  = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error   = $(job_name).$(jobid)
# Type de travail
# @ job_type= parallel
# Nombre de processus demandes
# @ total_tasks = 1024
# Temps ELAPSED max. pour l'ensemble du job en hh:mm:ss
# @ wall_clock_limit = 20:00:00
# @ queue
cd $LOADL_STEP_INITDIR
\rm -rf core*
\rm -rf *rsl.*
\rm -rf wrfi*
\rm -rf wrfo*
poe ./ideal.exe
mv rsl.error.0000 ideal.rsl.error.0000
mv rsl.out.0000 ideal.rsl.out.0000
poe ./wrf.exe
