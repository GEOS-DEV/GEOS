#!/bin/bash -l
#SBATCH --job-name=weibull2
#SBATCH --partition=bigmem
#SBATCH -A rhurley6_bigmem
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --export=ALL
#SBATCH --mail-type=end
#SBATCH --mail-user=sghosh29@jhu.edu

start=$(date +%s)

source /home/sghosh29/data-rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/llnl-env/bin/activate

cd /home/sghosh29/data_rhurley6/sohanjit/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/

./runClean_sghosh29.sh 1

#cd /data1/sghosh29/Working_MPM_LLNL/testGEOS/thin/

#mpirun -np 64 /data1/sghosh29/Working_MPM_LLNL/GEOS/build-system76-pc-clang-release/bin/geosx -i mpm_thin.xml -x 8 -y 8 -z 1

#python /data1/sghosh29/Working_MPM_LLNL/testGEOS/PostProcessing/thin.py

#cd /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull4/

#mpirun -np 48 /home/sghosh29/data_rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx -i /home/sghosh29/data_rhurley6/sohanjit/finalGEOS/weibull4/mpm_weibull4.xml -x 8 -y 6 -z 1

end=$(date +%s)
runtime=$((end-start))
echo "Total runtime: $runtime seconds"
