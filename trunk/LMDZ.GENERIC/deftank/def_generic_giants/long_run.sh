###Script to launch several successive jobs on Gnome 

# In your directory, you need :
# - this script
# - the reference script for Load Leveler, containing the settings of a job (which .def to use, which executable, Load Leveler settings...).
 

#!/bin/bash

##### MODIFY THIS :  
n=3  #total number of jobs 
ref='launch_128_96_64_ref' #name of the reference script 


#Following commands are specific to Load Leveler : change them if you use another batch scheduling system
ini_dir='$LOADL_STEP_INITDIR'		# Absolute path to the submission directory, as named in the reference script 
run_dir='$LOADL_STEP_INITDIR/${tmpdir}'	#Generic absolute path to the run directory, as named in the reference script (eg. job_name.job_id etc...)

sub_cmd=llsubmit 			#command to submit a job

##############################################################################################################################

# THE REFERENCE SCRIPT IS READ :

#we check that the lines about the copy of the startfi.nc follows the one with
#start.nc 

ligne="$(grep -n 'cp .*start.nc' ${ref} | cut -d: -f1)"
lignefi="$(grep -n 'cp .*startfi.nc' ${ref} | cut -d: -f1)"

if [ ${lignefi} -ne "$(expr ${ligne} \+ 1)" ]
then
	echo 'Error :'
	echo 'In '${ref}' ,the line about startfi.nc must follow the line about start.nc'
	exit
fi

#CREATION OF THE 1ST SCRIPT
cp ${ref} launch1

#CREATION OF THE LAST SCRIPT

#same as the reference script except that its start*.nc files are the restart*.nc files produced by the second to last simulation
 
head -n $(expr ${ligne} \- 1 ) launch1 > launch${n}
echo 'cp -rf ' ${ini_dir}'/zeLAST/restart.nc ' ${run_dir}'/start.nc' >> launch${n} 
echo 'cp -rf ' ${ini_dir}'/zeLAST/restartfi.nc ' ${run_dir}'/startfi.nc' >> launch${n}
tail -n +$(expr ${lignefi} \+ 1 ) launch1 >> launch${n}


#1ST SCRIPT IS COMPLETED AND THE OTHERS (from 2 to n-1 ARE CREATED) 

# Scripts 2 to n-1 are created from the last script (as jobs 2 to n-1 must start from restart*.nc produced by the previous job).
# Commands are added to Scripts 1 to n-1 to launch the next job.

for ((i=1 ; i<=${n}; i+=1))
do 	
	if [ ${i} -ne 1 ] && [ ${i} -ne ${n} ]
	then  
		cp launch${n} launch${i}
	fi 

	#Link toward the repertory of the last simulation
	echo 'cd ' ${ini_dir} >> launch${i}
	echo 'rm -f zeLAST' >> launch${i}
	echo 'lastdir=`ls -ltrd */ | awk' "'{print" '$NF}'"'" '| tail -n 1 `' >> launch${i}	
	echo 'ln -sf ${lastdir} zeLAST' >> launch${i}
	echo >> launch${i}
	
	if [ ${i} -ne ${n} ]
	then 
		#Launch the next job
		echo '#Lancement du job script suivant' >> launch${i}
		echo 'cd' ${ini_dir} >> launch${i}
		echo ${sub_cmd}' launch'$(expr $i \+ 1) >> launch${i}
	fi
done

 
#LAST RUN : CREATION OF outputs.txt
# concatnc.e < outputs.txt will concatenate all the diagfi.nc 

echo 'cd '${ini_dir} >> launch${n}
echo 'find . -maxdepth 2 -name "diagfi.nc" >> outputs.txt' >> launch${n}

echo 'echo >> '${ini_dir}'/outputs.txt'  >> launch${n}
echo 'echo 0 >> '${ini_dir}'/outputs.txt' >> launch${n}
echo 'echo ls >> '${ini_dir}'/outputs.txt' >> launch${n}
echo 'echo all >> '${ini_dir}'/outputs.txt'  >> launch${n}
echo 'echo >> '${ini_dir}'/outputs.txt'  >> launch${n}
 
#LAUNCH THE 1ST JOB
${sub_cmd} launch1

