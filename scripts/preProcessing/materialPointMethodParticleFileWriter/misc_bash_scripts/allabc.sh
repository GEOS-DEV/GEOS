#!/bin/bash -l
#SBATCH --job-name=allabc
#SBATCH --partition=parallel
#SBATCH -A rhurley6
#SBATCH --time=07-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --qos=extended
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

start=$(date +%s)

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/abc2/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_abc2.xml -x 16 -y 12 -z 1

cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/abc3/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_abc3.xml -x 16 -y 12 -z 1


cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/abc4/

mpirun -np 192 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i mpm_abc4.xml -x 16 -y 12 -z 1


end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
