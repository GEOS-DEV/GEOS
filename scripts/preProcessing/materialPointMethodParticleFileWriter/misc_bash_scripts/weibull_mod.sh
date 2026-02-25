#!/bin/bash -l
#SBATCH --job-name=weibull1
#SBATCH --partition=parallel
#SBATCH -A rhurley6
#SBATCH --time=03-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

# start=$(date +%s)

# source /home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/llnl-env/bin/activate

# cd /home/sghosh29/data_rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/

# fileNames=(
# 	#elasticBlockUni
# 	#planestrain
# 	#ceramicDamage
# 	#full3Dsmall
# 	#full3Dpart2
# 	#mesh4
# 	#normal4
# 	#abc4
# 	weibull1

# )
# # ==========================================================================================================================================
# # This should be the location of the input file and anything else you need to copy over:
# fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/weibull/'

# # This is where you want to run the simulation, which should be on a large parallel file system (lustre or workspace).  
# # This directory should exist.  sub-directory with fileName will be created
# runLocation='/home/sghosh29/data-rhurley6/sohanjit/finalGEOS/' # Quartz and Ruby

# # This gets the username so the right localdefs file can be used.
# userName="$(whoami)"

# # You shouldn't need to modify any of this.  It checks if the fileName you provided isn't a null string (so you don't delete the parent directory!)
# # then it copies files needed to run, changes directory, and runs the particle file writer.  That script will create a batch file and submit the job
# # if those options are enabled:
# #
# # For certain input options (like using CT data) you'll also need to copy over any other needed files (CT data, input tables, etc.)
# for fileName in "${fileNames[@]}"
# do
# 	if [ -n "$fileName" ]
# 	then
# 		if [ $# -eq 0 ] || [ $1 -eq 1 ];
# 		then
# 			num_tasks="1"
# 			echo "Running job on 1 process: ${fileName}"
# 		else
# 			num_tasks="$1"
# 			echo "Running job on ${num_tasks} processes: ${fileName}"
# 		fi 

# 		aborted=false
# 		if [ -d $runLocation/$fileName/ ] && [ -z "$SLURM_JOBID" ];
# 		then
# 			echo "Directory ${runLocation}/${fileName} exists."
# 			while true; do
# 				read -p "Do you wish to overwite? " yn
# 				case $yn in
# 					[Yy]* ) echo "Overwriting..."; break;;
# 					[Nn]* ) echo "Aborted overwrite..."; aborted=true; break;;
# 					* ) echo "Please answer yes (Y/y) or no (N/n).";;
# 				esac
# 			done
# 		fi

# 		if [ $aborted = true ]
# 		then
# 			continue
# 		fi

# 		rm -rf $runLocation/$fileName/											# delete old results for the same fileName!!!
# 		mkdir -p $runLocation/$fileName/										# create the run/output directory

# 		cp $fileLocation/pfw_input_$fileName.py $runLocation/$fileName          # copy the input file
# 		cp particleFileWriter.py $runLocation/$fileName           				# copy the preprocessor
# 		cp pfw_check.py $runLocation/$fileName                    				# copy the autoRestart script
# 		cp pfw_geometryObjects.py $runLocation/$fileName          				# copy the geometry object functions
# 		cp userDefs_$userName.py $runLocation/$fileName           				# copy the local path information	

# 		cd $runLocation/$fileName                                 				# move to the run location
# 		if [ $# -eq 0 ] || [ $1 -eq 1 ];
# 		then
# 			python3 particleFileWriter.py pfw_input_$fileName
# 		else
# 			mpirun -n ${num_tasks} python3 particleFileWriter.py pfw_input_$fileName          # launch the VML
# 		fi
# 		echo # Print empty line for legibility
# 	fi
# done

#cd /data1/sghosh29/Working_MPM_LLNL/testGEOS/thin/

#mpirun -np 64 /data1/sghosh29/Working_MPM_LLNL/GEOS/build-system76-pc-clang-release/bin/geosx -i mpm_thin.xml -x 8 -y 8 -z 1

#python /data1/sghosh29/Working_MPM_LLNL/testGEOS/PostProcessing/thin.py

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull1/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull1/mpm_weibull1.xml -x 16 -y 12 -z 1

rm -r mpm_weibull1_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r weibull1/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/


end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
