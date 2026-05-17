#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A 
#SBATCH --export=NONE
#SBATCH -p pdebug

echo "Launching srun command..."
export OMP_NUM_THREADS=1
srun --export=ALL,MV2_SMP_USE_CMA=0,PSM2_KASSIST_MODE=none -n 24 /usr/WS1/homel1/GEOS/scripts/preProcessing/materialPointMethodParticleFileWriter/examples/benchy/_suite_runs/benchy/benchy/geosx -i mpm_benchy.xml -x 4 -y 2 -z 3 
echo "srun command has completed, good bye."
