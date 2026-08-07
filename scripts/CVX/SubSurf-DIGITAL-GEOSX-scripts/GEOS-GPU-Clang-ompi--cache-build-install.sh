#!/bin/bash 
#
# Chevron generic clone-config-build-install script for GEOS and GEOSX on x86_64 platforms
#
# Compilers : GNU C/C++ Fortran and Nvidia NVCC (CUDA)
# MPI Stack : OpenMPI HPC-X 
# ARCH :      x86_64
# OS :        RHEL7/8
#
# Michael E. Thomadakis, Innovation and HPC R&D, CTC 
#
# NOTE :
# All references below to GEOSX should be intepreted as referring to GEOS as well
#
# 
# Setup Python, git, git-lfs, compilers
#
# Download and build the thirdPartyLibs for GEOS via git and git-lfs from
#   https://github.com/GEOS-DEV/thirdPartyLibs.git
#
# Download and build GEOS via git and git-lfs from
#   https://github.com/GEOS-DEV/GEOS.git
# Documentation
#  https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/QuickStart.html Steps



# Usage  
#
#  [sroot=src_dir][iroot=install_dir] [Pre ENVVARS] $0 CONFIG_NAME|--help 
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

: ${gitTPL:="https://github.com/GEOS-DEV/thirdPartyLibs.git"}
: ${gitGEOS:="https://github.com/GEOS-DEV/GEOS.git"}

# Old URL
#  https://github.com/GEOS/GEOS.git 
#  https://github.com/GEOSX/GEOSX.git 

export gitTPL gitGEOS 

: ${CUDA:="12"} 
: ${mpi:="hpcx"}

if [[ "$CUDA" == "12" ]]; then
    : ${gpu_comp:="CUDA/12.6"}
    : ${comp:="gcc/13.2.0-rh8"} 
else
    : ${gpu_comp:="CUDA/11.8"}
    : ${comp:="gcc/11.4.0-rh8"} 
fi

  ${comp}+=" llvm/20.1.6-gcc-14.2.0-ke2wmoz "

: ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf git/2.39.1 git/lfs_3.2.0 CMake_3.28.3 "}
#: ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf git/2.39.1 git/lfs_3.2.0 CMake_4.0.3 "}
#: ${modules:="bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf git/2.39.1 git/lfs_3.2.0 CMake_3.28.3 "}

# : ${modules:=" bison/3.8.2-gcc-13.2.0-io5a5qf flex/2.6.3-gcc-13.2.0-io5a5qf /devl/geophys/util/modules/ModuleFiles/git/2.27.0  CMake_3.24.1  "}
# : ${modules:=" /devl/geophys/util/modules/ModuleFiles/git/2.27.0  CMake_3.24.1 CUDA/UL_12.0 10.2.0 "}
# : ${modules:=" /devl/geophys/util/modules/ModuleFiles/git/2.27.0 CUDA/NVHPC_23.1.0 CMake_3.24.1 10.2.0 "}
# : ${modules:=" /devl/geophys/util/modules/ModuleFiles/git/2.27.0 CUDA/NVHPC_22.7.0 CMake_3.24.1 10.2.0 "}
# x86_64-12.2.0
# CMake_3.18.4  

export MODULES=" ${modules} ${mpi} ${comp} ${gpu_comp} "

# Check if we were sumitted as a PBS job
if [ -n "${PBS_O_WORKDIR}" ]; then
    cd ${PBS_O_WORKDIR}
    export PBS=1
    : ${scrdir:="$(cd . ; pwd)"}
    currdir=$(basename ${scrdir})
    if [ "$currdir" != "SubSurf-DIGITAL-GEOSX-scripts" ]; then 
	echo "## $0 Error : you may only submit a PBS batch job to config-build GEOSX scripts from within the SubSurf-DIGITAL-GEOSX-scripts"
	echo "## pwd(=$currdir) != SubSurf-DIGITAL-GEOSX-scripts; Exiting ..."
    	exit 10; 
    fi
    
    echo "## Running under PBS ; PBS_ENVIRONMENT=$PBS_ENVIRONMENT"
    echo "## switched back to ${PBS_O_WORKDIR} "
else
    scrdir="$(cd $(dirname $0) ; pwd)"
    currdir=$(basename ${scrdir})
    export PBS=0
fi

export scrdir currdir

# Default location of Chevron's host-config .cmake 
: ${HOST_CONFIG_DIR:="${scrdir}/GEOS/host-configs/CVX"}
export HOST_CONFIG_DIR

if [ $# -lt 1 ]; then 
    export general_help=1 
    source ${scrdir}/GEOS-quick-help.sh
    echo "## Existing CVX CONFIG_NAME.cmake at ${HOST_CONFIG_DIR}; set CVX_CONFIG_DIR to other location if desired "
    for hcf in $(ls -1 ${HOST_CONFIG_DIR}/*GPU*ompi*.cmake ${HOST_CONFIG_DIR}/*GPU*MPI*.cmake); do 
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
export GEOS_GPU_WORKAROUND=1 

#  Initialize CVX Python
unalias python python3 >& /dev/null
export ANACONDA_VER="2024.02"
# export ANACONDA_VER="2023.09"

source /util/Anaconda/${ANACONDA_VER}/etc/profile.d/conda.sh
conda activate 
export CVX_PYTHON3=/util/Anaconda/${ANACONDA_VER}/bin/python
# export PYTHON3=$(which python3) 

# Initialize Chevron MODULES system
source /data/saet/hpcrnd/utils/bin/modules_init.sh 
# WHO=${USER} source /data/saet/hpcrnd/utils/bin/modules_init.sh 

# Make git-lfs visible via PATH
# GIT_LFS="/data/saet/mtml/software/x86_64/bin"
# export PATH="${GIT_LFS}:${PATH}"
 
module load ${MODULES}
module list

# Determine host compiler and version
export GCC_VER=$(gcc --version | gawk '$1=="gcc" && $2=="(GCC)" {print $3 }')
export GCCVER=$GCC_VER
export HOST_COMPILER_VER="$(clang++ --version | gawk '$1=="clang" && $2=="version" {print $3 }')"
export HOST_COMPILER="CLANG_${HOST_COMPILER_VER}"

export GCC_TOOLCHAIN=$(dirname $(dirname $(which gcc)))

# Set compilers for mpi wrappers 
export OMPI_CC=$(which clang)
export OMPI_CXX=$(which clang++)
export OMPI_F77=$(which gfortran)
export OMPI_FC=$(which gfortran)

# Compilers 
export GEOSX_CC=${OMPI_CC}
export GEOSX_CXX=${OMPI_CXX}
export GEOSX_FORT=${OMPI_FC}

# Open MPI specifics
export GEOSX_MPI_DIR=$HPCX_MPI_DIR

export GEOSX_MPIRUN=${GEOSX_MPI_DIR}/bin/mpirun
export GEOSX_MPICC=${GEOSX_MPI_DIR}/bin/mpicc
export GEOSX_MPICXX=${GEOSX_MPI_DIR}/bin/mpic++
export GEOSX_MPIFORT=${GEOSX_MPI_DIR}/bin/mpifort

export MPIRUN=$GEOSX_MPIRUN

# MPI independent logic
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
echo "## scrdir= $scrdir ;  currdir=$currdir "

echo "# Chevron HPC Innovation and R&D generic clone-cache-config-build-install script for GEOS / TPL
#
# Compilers : GCC $GCC_VER (modulefile= ${comp}, gpu_comp= $gpu_comp ; CUDA= $CUDA)
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
