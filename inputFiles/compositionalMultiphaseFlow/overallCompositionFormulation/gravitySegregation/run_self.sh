source /util/hpcx/hpcx-v2.14-gcc-MLNX_OFED_LINUX-5-redhat8-cuda11-gdrcopy2-nccl2.16-x86_64/hpcx-init.sh
source /data/saet/hpcrnd/utils/bin/modules_init.sh
hpcx_load
which mpirun
 
export GEOSX_DIR=/data/rpo_fjvz/codes/GEOS/build-self_build-relwithdebinfo

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GEOSX_DIR/lib
 
# Show loaded modules and location of geosx binary
module list
 
# Show current directory
pwd
 
export OMP_NUM_THREADS=1

mpirun -np 1 $GEOSX_DIR/bin/geosx -i ./grav_seg_drain_Z.xml -o ./resultsZ > grav_seg_drain_Z.out
hpcx_unload
