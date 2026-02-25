#!/bin/bash -l
#SBATCH --job-name=3DRBC
#SBATCH --partition=bigmem
#SBATCH -A rhurley6_bigmem
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --qos=extended
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu


start=$(date +%s)

./runClean_sghosh29.sh 1

echo "Getting Strength Distribution"

python /home/sghosh29/data_rhurley6/sohanjit/testGEOS/PostProcessing/3Ddist.py

end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"


#echo "Launching srun command..."
#srun -n 48 /home/sghosh29/data-rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i /home/sghosh29/data_rhurley6/sohanjit/testGEOS/2D_HighProjectileY/mpm_2D_HighY.xml -x 8 -y 6 -z 1
#echo "srun command has completed, good bye."

