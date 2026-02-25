#!/bin/bash -l
#SBATCH --job-name=3Dsmall_SS
#SBATCH --partition=bigmem
#SBATCH -A rhurley6_bigmem
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --qos=extended
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

# start=$(date +%s)

# source /home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/llnl-env/bin/activate

# cd /home/sghosh29/data_rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter

# fileNames=(
# 	3Dsmall
	

# )
# # ==========================================================================================================================================
# # This should be the location of the input file and anything else you need to copy over:
# #fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/Ftable/'
# #fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/planestrain/'
# #fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/ceramicDamage/'
# #fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/full3D_small/'
# #fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/full3D_part2/'
# fileLocation='/home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/3Dsmall/'

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



cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/3Dsmall/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/3Dsmall/mpm_3Dsmall.xml -x 8 -y 6 -z 4

rm -r mpm_3Dsmall_restart*

sshpass -p "iittojhu2021!" scp -r vtkOutput* sghosh29@10.161.161.71:/data3/sohanjit/3DGEOS/3Dsmall_revised/

end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
