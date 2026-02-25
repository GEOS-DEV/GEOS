#!/bin/bash -l
#SBATCH --job-name=allnormal
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

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/normal/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_normal.xml -x 16 -y 12 -z 1
rm -r mpm_normal_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r normal/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/normal1/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_normal1.xml -x 16 -y 12 -z 1
rm -r mpm_normal1_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r normal1/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/normal2/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_normal2.xml -x 16 -y 12 -z 1
rm -r mpm_normal2_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r normal2/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/normal3/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_normal3.xml -x 16 -y 12 -z 1
rm -r mpm_normal3_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r normal3/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/normal4/
mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_normal4.xml -x 16 -y 12 -z 1
rm -r mpm_normal4_restart_*
cd ..
sshpass -p "iittojhu2021!" scp -r normal4/ sghosh29@10.161.161.71:/data1/sghosh29/Working_MPM_LLNL/finalGEOS/

end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
