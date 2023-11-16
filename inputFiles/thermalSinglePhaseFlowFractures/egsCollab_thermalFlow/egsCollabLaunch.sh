#!/bin/tcsh
##### These lines are for Slurm
#SBATCH -N 2
#SBATCH -J egsCollab
#SBATCH -t 10:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A vortex
#SBATCH -o /usr/workspace/cusini1/geosx/geosx_integratedTestsRebase/GEOSX/egsCollab.out

##### These are shell commands
date
cd /usr/workspace/cusini1/geosx/geosx_integratedTestsRebase/GEOSX/build-quartz-gcc-12-release
##### Launch parallel job using srun
setenv OMP_NUM_THREADS 1
srun -N 2 -n64 ./bin/geosx -i ../inputFiles/thermalSinglePhaseFlowFractures/egsCollab_thermalFlow/egsCollab_thermalFlow_injection_fine.xml -x 4 -y 4 -z 4 -o Output_collab
echo 'Done'