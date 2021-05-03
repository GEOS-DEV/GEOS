#!/bin/sh
SBATCH --job-name=test_pygeosx
#SBATCH --time=00:15:00
SBATCH -n 8
SBATCH -N 1
SBATCH -C bora
#SBATCH --exclusive
#echo ${SLURM_JOB_ID} > IDseg_large8Hz${SLURM_JOB_ID}.txt
#RUNDIR=/beegfs/mndiaye/GEOSXDATA_Benchmarks/676x676x210_results/seg_large8Hz${SLURM_JOB_ID}

mkdir ${RUNDIR}
cd ${RUNDIR}
mkdir SISMO

mpirun </dev/null --map-by node python /beegfs/jbesset/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/pygeosx/main.py -i /beegfs/jbesset/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/benchmarks/pygeosx_test.xml -x 4 -y 2  
mv sismo* SISMO/
mv vtkOutpu* SISMO/
