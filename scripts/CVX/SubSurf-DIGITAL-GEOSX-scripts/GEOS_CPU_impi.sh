#!/bin/bash 
#
# Generic gbatch script to run installed GEOS on x86_64 platforms (any scheduler)
#
# Compilers : GCC 
# MPI Stack : IntelMPI OneAPI
# ARCH :      x86_64
#
# Michael E. Thomadakis, Innovation and HPC R&D, Chevron Technology Center
#
# 
# Sets up enviroment modules, Python, git, git-lfs, compilers and the GEOSX variant 
#


# Usage  
# NP=np PPN=ppn [OMP_NUM_THREADS=ont] {INITRC=/full/path/init/XX-impiXX.rc|GESX_DIR=/full/pathto/GEOS} $0 GEOS arg list ..... 

if [ $# -lt 1 ]; then
    echo "## Usage : NP=np PPN=ppn OMP_NUM_THREADS=ont [INITRC=initrc|GESX_DIR=/full/pathto/GEOS] $0 GEOS arg list .....  "
    exit -10;
fi

[ -n "$DEBUG" ] && set -x 

: ${V:="0"}
export V

: ${NP:="1"}
export NP

: ${PPN:="${NP}"}
export PPN

export H=$(hostname -s)

export EXPD="."

export B0=$(basename $0)

## Prepare the directories and hostfiles 
# Check if we were sumitted as a PBS job
if [ -n "$SLURM_NODELIST" ]; then
    scrdir="$(cd "." ; pwd)"
    # scrdir="$(cd $(dirname $0) ; pwd)"
    export PBS=0
    expname="slurm-${SLURM_JOBID}";
    HOSTFILE="${EXPD}/${expname}.hosts";
    scontrol show  hostnames $SLURM_NODELIST > ${HOSTFILE}
    echo "## $B0 : Running under SLURM (${SLURM_JOBID}) at $(pwd)"
elif [ -n "${PBS_O_WORKDIR}" ]; then
    cd ${PBS_O_WORKDIR}
    scrdir="$(cd "." ; pwd)"
    export PBS=1
    HOSTFILE=${PBS_NODEFILE} ;
    expname="$(/bin/basename ${PBS_NODEFILE})"
    printf "## $B0 : Running under PBS ($PBS_JOBID); switched now at ${PBS_O_WORKDIR}\n## (pwd=$(pwd)) \n"
else 
    echo " $0 Warning : This is not being run under a batch system. Assumes a local host run."
    scrdir="$(cd $(dirname $0) ; pwd)"
    HOSTFILE="${EXPD}/${H}.hosts" 
    echo $H > $HOSTFILE
    expname=$(/bin/basename $HOSTFILE ".hosts");
    echo "## $B0 : Not being run under a batch system at $(pwd) "
fi

currdir=$(basename $(pwd))

export expname scrdir currdir
export HOSTFILE

HOSTTAB="${EXPD}/${expname}.hosttab"
gawk '$1 !~ /#/ && NF > 0' $HOSTFILE  | sort -u > $HOSTTAB
export N_nodes=$(wc $HOSTTAB | gawk '{print $1}'); 
export HOSTFILE HOSTTAB N_nodes
echo "# HOSTFILE= $HOSTFILE HOSTTAB= $HOSTTAB ; N_nodes= $N_nodes "
echo "## scrdir= $scrdir ;  currdir=$currdir "
export HOSTTAB

## Prepare the s/w stack 
: ${mpi:="mpi"}
: ${comp:="gcc/10.4.0"}
: ${modules:="/devl/geophys/util/modules/ModuleFiles/git/2.27.0 CMake_3.24.1 "}

export MODULES=" ${mpi} ${modules} ${comp} "

# Check for INITRC
if [ -n "$INITRC" ]; then 
    echo "## $0 Note : initializing with INITRC= $INITRC "
    echo "## Ignoring MODULES= $MODULES "
    source $INITRC
elif [ -n "$GEOSX_DIR" ]; then 
    export PATH="${GEOSX_DIR}/bin:${PATH}"
    export LD_LIBRARY_PATH="${GEOSX_DIR}/lib:${LD_LIBRARY_PATH}"
    # Initialize Chevron MODULES system
    source /data/saet/hpcrnd/utils/bin/modules_init.sh 
    # WHO=${USER} source /data/saet/hpcrnd/utils/bin/modules_init.sh 

    # Make git-lfs visible via PATH
    export PATH="/util/git-lfs/prod:${PATH}"
 
    echo "## Loading MODULES= $MODULES "
    module load ${MODULES}
    module list
else
    echo "## $0 Error : INITRC does not exist and GEOSX_DIR is underfined ! Exiting "
    exit 10; 
fi

export GEOSX_ARGS="$@"
export GEOSX=$(which geosx)
export GEOSX_MPIRUN=$(which mpirun)

: ${ipin_domain:="compact"}
: ${RANKMAP:="spread"} 

if [ -n "$OMP_NUM_THREADS" ]; then 
    export PE=$OMP_NUM_THREADS
    export MPIOPTS=" -genv I_MPI_PIN_DOMAIN ${PE}:${ipin_domain} -genv I_MPI_PIN_ORDER ${RANKMAP} "
else
    export PE=1
    export MPIOPTS=" -genv I_MPI_PIN_PROCESSOR_LIST allcores:map=${RANKMAP}"
fi
export MPIRUNOPTS=" -f $HOSTTAB -np $NP -ppn ${PPN} ${MPIOPTS} ${CLIMPIOPTS} "

if [ $V -gt 0 ]; then 
    echo "# \$0= $0, pwd= $(pwd)" 
    echo "# \$#= $#, \$*= [$*]" 
    echo "# INITRC= $INITRC "
    echo "# NP= $NP, PPN= $PPN, OMP_NUM_THREADS= $OMP_NUM_THREADS "
    echo "# MPIRUN= $MPIRUN "
    echo "# MPIRUNOPTS= $MPIRUNOPTS "
    echo "# geosx : $GEOSX "
    echo "# GEOSX_ARGS= [$GEOSX_ARGS] "
    echo "# Loaded modules : "
    module list >&2
fi

# Command line
if [ $V -gt 1 ]; then 
    which mpirun
    echo "# Batch Env Vars : "
    env | egrep 'SLURM|PBS|MPI|OMP|GEOS' | sort | gawk '{printf "# %s \n", $0}'
    lscpu 
fi

echo "# -------------------------------------------------------------------------
# Running GEOS on host=$(hostname) in=$(pwd) at $(date +%F-%H%M%S) as :
$PRE_MPIEXEC $GEOSX_MPIRUN $MPIRUNOPTS $GEOSX $GEOSX_ARGS"

$PRE_MPIEXEC $GEOSX_MPIRUN $MPIRUNOPTS $GEOSX $GEOSX_ARGS
echo "# -------------------------------------------------------------------------"
