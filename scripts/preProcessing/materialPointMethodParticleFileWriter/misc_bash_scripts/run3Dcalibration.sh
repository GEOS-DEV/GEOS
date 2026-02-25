#!/bin/bash -l
#SBATCH --job-name=3DRBC_0.19
#SBATCH --partition=bigmem
#SBATCH -A rhurley6_bigmem
#SBATCH --time=07-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --qos=extended
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

start=$(date +%s)

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/3DRBC/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/3DRBC/mpm_3DRBC.xml -x 8 -y 6 -z 4

rm -r mpm_3DRBC_restart*

sshpass -p "iittojhu2021!" scp -r vtkOutput2* sghosh29@10.161.161.71:/data3/sohanjit/3DGEOS/

end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
