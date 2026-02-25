#!/bin/bash -l
#SBATCH --job-name=allweibull
#SBATCH --partition=parallel
#SBATCH -A rhurley6
#SBATCH --time=03-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

start=$(date +%s)

# cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull1/
# mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_weibull1.xml -x 16 -y 12 -z 1
# rm -r mpm_weibull1_restart_*
# cd ..
# sshpass -p "iittojhu2021!" scp -r weibull1/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull2/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_weibull2.xml -x 16 -y 12 -z 1
rm -r mpm_weibull2_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r weibull2/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull3/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_weibull3.xml -x 16 -y 12 -z 1
rm -r mpm_weibull3_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r weibull3/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull4/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_weibull4.xml -x 16 -y 12 -z 1
rm -r mpm_weibull4_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r weibull4/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/


end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
