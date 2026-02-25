#!/bin/bash -l
#SBATCH --job-name=allrate
#SBATCH --partition=bigmem
#SBATCH -A rhurley6_bigmem
#SBATCH --time=02-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

start=$(date +%s)

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/rate1/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_rate1.xml -x 16 -y 12 -z 1
rm -r mpm_rate1_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r rate1/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/rate2/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_rate2.xml -x 16 -y 12 -z 1
rm -r mpm_rate2_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r rate2/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/rate3/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_rate3.xml -x 16 -y 12 -z 1
rm -r mpm_rate3_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r rate3/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/


end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
