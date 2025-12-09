#! /bin/bash
#
# CMake hypredrive builds
#
# MIchael E. Thomadakis HPC INnovation and R&D 
#

function echoc0 {
    echo "# $0 : $* "
}

function nowAt() {
    echo "# $1 Now at $(pwd), should be == $2 "
}

function cdnowAt() {
    cd $2
    echo "# $1 Now at $(pwd), should be == $2 "
}

function RCwarning() {
    if [[ $2 > 0 ]]; then
	echo "### $0 Warning : RC=$2; $3 "
    fi
}

function RCerror() {
    if [[ $2 > 0 ]]; then
	echo "### $0 ERROR ! RC=$2; $3; Exiting! "
	exit $RC;
    fi
}

: ${V:=0}

umask 0022

if [[ $V -ge 4 ]]; then 
    set -x
fi

: ${CLONE:="0"} ;
export CLONE

: ${HYPRE:="0"} ;
export HYPRE

: ${INIT_ONLY:="0"}
export INIT_ONLY 

: ${GPU_BUILD:="0"} 
export GPU_BUILD

SCRP=$(readlink -f $0)
SCRD=$(dirname $SCRP)
SCRB=$(basename $SCRP)
CWD=$(pwd)

# --- Parse options ---
: ${ENABLE_CUDA:=OFF}
: ${ENABLE_HIP:=OFF}
: ${UMPIRE:=0}
: ${CUDA_ARCH:="80"};
export CUDA_ARCH

: ${CMAKE_VERBOSE_MAKEFILE:="OFF"}
export CMAKE_VERBOSE_MAKEFILE

: ${BUILD_TYPE:="RelWithDebInfo"}
# BUILD_TYPE="Release"

: ${PKG:="hypredrive"}
: ${PKG_ver:="0.1.0"}

: ${ARCH:="x86_64"} ;
: ${OS:="RHEL8"} ;

: ${PUBLISH_INIT:="0"} ;
export PKG ARCH OS

if [[ $V -gt 0 ]]; then
    export VERBOSE=1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --relwithdebinfo)
      BUILD_TYPE="RelWithDebInfo"
      shift
      ;;
    --release)
      BUILD_TYPE="Release"
      shift
      ;;
    --debug)
      BUILD_TYPE="Debug"
      shift
      ;;
    --cuda)
      ENABLE_CUDA=ON
      SUFFIX="cuda"
      GPU_BUILD=1
      shift
      ;;
    ---hip)
      ENABLE_HIP=ON
      SUFFIX="hip"
      GPU_BUILD=1
      shift
      ;; 
    *) 
      echo "Unknown option: $1"
      echo "Usage: $0 [--cuda | --hip | --relwithdebinfo | --release | --debug ]"
      exit 1
      ;;
  esac
done

# Sanity Check
if [[ ${ENABLE_CUDA} == "ON" && ${ENABLE_HIP} == "ON" ]]; then
    echo -e "Please select either CUDA or HIP"
    exit 1
fi

build_type="$(echo ${BUILD_TYPE} | tr '[:upper:]' '[:lower:]')"
SUFFIX+="-${build_type}"

echoc0 "----------------------------------------------------------------------------"
echoc0 "SCRP=$SCRP"
echoc0 "SCRD=$SCRD"
echoc0 "SCRP=$SCRP"
echoc0 "CMAKE=$CMAKE CLONE=$CLONE HYPRE=$HYPRE UMPIRE=${UMPIRE}  "

if [ "${USER}" != "mtml" ]; then 
    if [ -z "$sroot" ]; then 
	echo "## Please set sroot='/file/system/dir for the location to clone sources and build codes' "
	export sroot=$(pwd) 
	echo "## setting sroot=${sroot}"
    fi
    if [ -z "$iroot" ]; then 
	echo "## Please set iroot='/file/system/dir for the location to install binaries and libraries'"
	export iroot=$(pwd)
	echo "## setting iroot=${iroot}"
    fi
else
    : ${sroot:="/data/saet/${USER}/src"} ; 
    : ${iroot:="/data/saet/hpcrnd/software"} ;
fi
: ${croot:="/dev/shm/${USER}/src"} ; 

if [[ ! -d $sroot ]]; then
    mkdir -p $sroot ; RC=$?
    if [ $RC -ne 0 ]; then
	echo "## $0 Warning : I cannot create $sroot ! Setting sroot = $(pwd) "  ;
	export sroot=$(pwd) 
    else
	echo "## $0 Note : Created sroot=$sroot "  ;
    fi
fi

sroot=$(readlink -f $sroot) ;
echo "## $0 Note : Using sroot=$sroot "  ;

 echo "## $0 Note : /data/saet/software installation "  ;
 if [[ ! -d $iroot ]]; then
     mkdir -p $iroot ; RC=$?
     if [ $RC -ne 0 ]; then
	 echo "## $0 Warning : I cannot create $iroot ! Setting iroot = $(pwd) "  ;
	 export iroot=$(pwd) 
     fi
 fi
echo "## $0 Note : Using iroot=$iroot "  ;
	
export sroot iroot croot
# 
: ${projdir:=""} 

# if [ "${projdir}" != "" ] ; then
#     export projdir="-${projdir}"
# fi

PKG="hypredrive"
PKG_ver="0.1.0"
PKG_git=" https://github.com/hypre-space/${PKG}.git"
export PKG PKG_ver PKG_git

if [[ ! -r $INITRC ]]; then
    echo "# $0 Error : INITRC ($INITRC) is undefined or initialization file does not exist  ! Exiting ... "
    exit 10;
fi

source $INITRC;
# HOST_CONFIG_FILENAME points to the host-configs/cmake file used
echoc0 "HOST_CONFIG_FILENAME=$HOST_CONFIG_FILENAME  "

module load automake/1.16.5-gcc-13.2.0-tibgshu autoconf/2.72-none-none-2enf465 libtool/2.4.7-gcc-13.2.0-xnegjiy 

echo "## $0 : Loaded modelefiles "
module list -v 

: ${ARCH:="x86_64"} ;
: ${OS:="RHEL8"} ;

export ARCH OS

cdnowAt $0 $sroot 

DT=$(date +%F_%H%M%S)

PKG_src="$(pwd)/${PKG}"
PKG_build="$PKG_src/${GEOSX_CONFIG_NAME}"
PKG_install="${iroot}/${ARCH}/${OS}/${PKG}/${PKG_ver}/${GEOSX_CONFIG_NAME}"

build_type="$(echo ${BUILD_TYPE} | tr '[:upper:]' '[:lower:]')"

STACK="${HOST_CONFIG_BASE}-${projdir}-${build_type}"
BUILDDIR="build-${STACK}"
INSTALLDIR="install-${STACK}"
SUFFIX=$INSTALLDIR
export build_type STACK BUILDDIR INSTALLDIR

# --- Umpire (optional dependency) ---
: ${CUDA_ARCH:="80"};
export CUDA_ARCH

UMPIRE_INSTALL=""
if [[ ( ${ENABLE_CUDA} == "ON" || ${ENABLE_HIP} == "ON" ) && $HYPRE == "1" ]]; then
    : ${UMPIRE:=1} 
    export UMPIRE
fi

export VERBOSE=1 V=1

    echoc0 "----------------------------------------------------------------------------"
    if [[ $GPU_BUILD == "1" ]]; then 
	hypredrive_build_cmake="${SCRD}/hypredrive_build_cmake_GPU.sh"
    else
	hypredrive_build_cmake="${SCRD}/hypredrive_build_cmake_CPU.sh"
    fi
    [[ -r ${hypredrive_build_cmake} ]] || RCerror $0 $? "Script ${hypredrive_build_cmake} to build via Cmake does not exist "
    echoc0 "Sourcing $hypredrive_build_cmake at $(pwd)" 
    source $hypredrive_build_cmake "$@" |& tee CMake-all_${DT}.log

echoc0 "-- END ---------------------------------------------------------------------"

if [ $V -ge 4 ]; then 
    set +x
fi
