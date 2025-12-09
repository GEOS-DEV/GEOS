#!/bin/bash 
#
# Chevron generic clone-config-build-install script for GEOS and GEOSX 
#
# Compilers : GCC 
# MPI Stack : IntelMPI 1-API or Classic
# OS :        RHEL7
#
#
# Michael E. Thomadakis, Innovation and HPC R&D, CTC 
#
#
# NOTE :
# All references below to GEOSX should be intepreted as referring to GEOS as well
#
# 
# Setup Python, git, git-lfs, compilers
#
# Download and build the thirdPartyLibs for GEOSX via git and git-lfs from
#   https://github.com/GEOSX/thirdPartyLibs.git
#
# Download and build GEOSX via git and git-lfs from
#   https://github.com/GEOSX/GEOSX.git
# Documentation
#  https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/QuickStart.html Steps
#PBS -S /bin/bash
 
# Usage  
#
#   [sroot=src_dir][iroot=install_dir] [Pre ENVVARS] $0 CONFIG_NAME|--help  
#
# Set CONFIG_NAME to a matching conf.cmake in GEOSX/host-configs/CVX
# 

#
## ------------------------------------------------------------------------------------
## Check and prepare for a batch run
## ------------------------------------------------------------------------------------

: ${WHO:=${USER}}

: ${V:="0"} ; 
export V

export H=$(hostname)

: ${BUILDID:=$(date +%F-%H%M%S)}
export BUILD_DT=$BUILDID 

: ${gitTPL:="https://github.com/GEOS-DEV/thirdPartyLibs.git"}
: ${gitGEOS:="https://github.com/GEOS-DEV/GEOS.git"}

# Old URL
#  https://github.com/GEOS/GEOSX.git 
#  https://github.com/GEOSX/GEOSX.git 

export gitTPL gitGEOS 

: ${CUDA:="12"} 

: ${mpi:="intel/oneapi/mpi/latest"}
# : ${mpi:="mpi/2021.5.1 UCX_1.13.0"}
# : ${mpi:="mpi UCX_1.13.0"}
: ${comp:="gcc/14.2.0-rh8"}
# : ${comp:="gcc/13.2.0-rh8"}
: ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf git/2.39.1 git/lfs_3.2.0 CMake_3.28.3 "}
#: ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf git/2.39.1 git/lfs_3.2.0 CMake_4.0.3 "}
# : ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf /devl/geophys/util/modules/ModuleFiles/git/2.27.0 git-lfs_3.2.0 CMake_3.28.3 "}

export MODULES=" ${mpi} ${modules} ${comp} "

# Check if we were sumitted as a PBS job
if [ -n "${PBS_O_WORKDIR}" ]; then
    cd ${PBS_O_WORKDIR}
    export PBS=1
    : ${scrdir:="$(cd . ; pwd)"}
    currdir=$(basename ${scrdir})
    if [ "$currdir" != "SubSurf-DIGITAL-GEOSX-scripts" ]; then 
	echo "## $0 Error : you may submit a GEOSX config-build script as a PBS batch job ONLY from within the SubSurf-DIGITAL-GEOSX-scripts"
	echo "## pwd(=$currdir) != SubSurf-DIGITAL-GEOSX-scripts; Exiting ..."
    	exit 10; 
    fi
    echo "## Running under PBS ; PBS_ENVIRONMENT=$PBS_ENVIRONMENT"
    echo "## Switched back to ${PBS_O_WORKDIR} "
else
    scrdir="$(cd $(dirname $0) ; pwd)"
    export PBS=0
    currdir=$(basename ${scrdir})
fi

export scrdir currdir

# Default location of Chevron's host-config .cmake 
: ${HOST_CONFIG_DIR:="${scrdir}/GEOS/host-configs/CVX"}
export HOST_CONFIG_DIR

if [ $# -lt 1 ]; then 
    export general_help=1 
    source ${scrdir}/GEOS-quick-help.sh
    echo "## Existing CVX CONFIG_NAME.cmake at ${HOST_CONFIG_DIR}; set CVX_CONFIG_DIR to other location if desired "
    for hcf in $(ls -1 ${HOST_CONFIG_DIR}/*CPU*impi*.cmake); do 
	hc=$(basename $hcf .cmake)
	printf " %54s  (%s)\n" $hc $hcf 
    done

    exit 10; 
fi

export CONFIG_NAME=$1

if [ "${CONFIG_NAME}" == "--help" ]; then 
    export general_help=1 
    source ${scrdir}/GEOS-quick-help.sh

    less ${scrdir}/README.md
    exit 1;
fi

hcn="${HOST_CONFIG_DIR}/$(basename $CONFIG_NAME .cmake).cmake"

if [ ! -r ${hcn} ]; then
    echo "## $0 Error : CONFIG_NAME file $hcn does note exist! "
    exit 2; 
fi


# -------------------------------------------------------------------------- #
## Initialize BULD environment 
# export GEOS_GPU_WORKAROUND=0 

#  Initialize CVX Python
unalias python python3 >& /dev/null
export ANACONDA_VER="2024.02"
# export ANACONDA_VER="2023.09"

source /util/Anaconda/${ANACONDA_VER}/etc/profile.d/conda.sh
conda activate 

export CVX_PYTHON3=/util/Anaconda/${ANACONDA_VER}/bin/python

# Initialize Chevron MODULES system
source /data/saet/hpcrnd/utils/bin/modules_init.sh 
# WHO=${USER} source /data/saet/hpcrnd/utils/bin/modules_init.sh 

# Make git-lfs visible via PATH
# export PATH="/util/git-lfs/prod:${PATH}"
 
module load ${MODULES}
module list

# Determine host compiler and version
export GCC_VER=$(gcc --version | gawk '$1=="gcc" && $2=="(GCC)" {print $3 }')
export GCCVER=$GCC_VER
export HOST_COMPILER_VER="$(clang++ --version | gawk '$1=="clang" && $2=="version" {print $3 }')"
export HOST_COMPILER="CLANG_${HOST_COMPILER_VER}"

# Set compilers for mpi wrappers 
export I_MPI_CC=$(which gcc)
export I_MPI_CXX=$(which g++)
export I_MPI_FC=$(which gfortran)
export I_MPI_F90=$(which gfortran)

# Compilers 
export GEOSX_CC=${I_MPI_CC}
export GEOSX_CXX=${I_MPI_CXX}
export GEOSX_FORT=${I_MPI_FC}

# MPI independent logic
export GEOSX_MPI_DIR=$I_MPI_ROOT

if [[ $I_MPI_ROOT =~ "mpi/202" ]]; then 
    export GEOSX_MPIRUN=${GEOSX_MPI_DIR}/bin/mpirun
    export GEOSX_MPICC=${GEOSX_MPI_DIR}/bin/mpiicc
    export GEOSX_MPICXX=${GEOSX_MPI_DIR}/bin/mpiicpc
    export GEOSX_MPIFORT=${GEOSX_MPI_DIR}/bin/mpiifort
else
    export GEOSX_MPIRUN=${GEOSX_MPI_DIR}/intel64/bin/mpirun
    export GEOSX_MPICC=${GEOSX_MPI_DIR}/intel64/bin/mpiicc
    export GEOSX_MPICXX=${GEOSX_MPI_DIR}/intel64/bin/mpiicpc
    export GEOSX_MPIFORT=${GEOSX_MPI_DIR}/intel64/bin/mpiifort
fi

export MPIRUN=$GEOSX_MPIRUN

export MPI=$($MPIRUN -version | gawk '
      { if ($1=="Intel(R)" && $2=="MPI") {printf "impi"}
        else if ($2=="(Open" && $3=="MPI)") {printf "ompi"} }' )

if [ "${MPI}" == "impi" ]; then
    export MPI_VER=$($MPIRUN -version | gawk '$1=="Intel(R)" && $2=="MPI" {printf "impi-%d.%d",$8, $10}')
else
    export MPI_VER=$($MPIRUN -version | gawk '$2=="(Open" && $3=="MPI)" {printf "ompi-%s",$4}')
fi

echo "## $0 : Running on $H; pwd=$(pwd) "
echo "# CONFIGNAME=${CONFIG_NAME} "
echo "## Scripts dir scrdir= $scrdir ;  currdir=$currdir "

echo "## Chevron HPC Innovation and R&D generic clone-cache-config-build-install script for GEOS 
#
# Compilers : GCC $GCC_VER (modulefile= ${comp})
# MPI Stack : $MPI $MPI_VER (modulefile= ${mpi})
"

# -------------------------------------------------------------------------- #
# Beyond this point IntelMPI and OpenMPI should be identical
# -------------------------------------------------------------------------- #

BESCRIPT="$scrdir/GEOS--cache-build-install.sh"

if [ ! -r $BESCRIPT ]; then
    echo "## $0 Error : Backend script $BASCRIPT does not exist! Exiting ... "  ;
    exit 30
fi

source $BESCRIPT $CONFIG_NAME

# -------------------------------------------------------------------------- #
