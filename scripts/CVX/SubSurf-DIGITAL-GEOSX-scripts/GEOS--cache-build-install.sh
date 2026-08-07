#!/bin/bash 
#
# Chevron generic clone-config-build-install framework for GEOS and GEOSX 
#
# Backend that 
#   1. preps user environemt
#   2. clones-configures-builds and installs 
# thirPartyLibs and GEOS
#
#
# Michael E. Thomadakis, Innovation and HPC R&D, CTC 
#
# NOTE :
# All references below to GEOSX should be intepreted as referring to GEOS as well
#
#
# Download and build the thirdPartyLibs for GEOS via git and git-lfs from
#   https://github.com/GEOS-DEV/thirdPartyLibs.git
#
# Download and build GEOS via git and git-lfs from
#   https://github.com/GEOS-DEV/GEOS.git
# Documentation
#  https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/QuickStart.html Steps


# Usage $0 CONFIG_NAME  
#
#  sourced fron script in $scrdir to carry out the work 
#
# Set CONFIG_NAME to a matching conf.cmake in GEOS/host-configs/CVX
# 

#
## ------------------------------------------------------------------------------------

# -------------------------------------------------------------------------- #
# Beyond this point IntelMPI and OpenMPI should be identical
# -------------------------------------------------------------------------- #

function RCwarning() {
    if [[ $2 > 0 ]]; then
	echo "### $0 Warning : RC=$2; $3 "
    fi
}

function RCerror() {
    if [[ $2 > 0 ]]; then
	echo "### $0 ERROR ! RC=$2; Exiting; $3 "
	exit $RC;
    fi
}

echo "## $0 Note : running backend script at $(pwd) for CONFIG_NAME=$CONFIG_NAME " ;

umask 0022

: ${V:="0"} ; 
export V

if [ $V -gt 1 ]; then
    export VERBOSE=1
fi

[ "$DEBUG" == "1" -o $V -ge 10 ] && set -x

: ${GEOS_VER:="1.1.0"} ; 
export GEOS_VER

: ${CHAP:="0"} ; 
export CHAP

: ${N_cores:=16}; 
: ${NUM_PROC:=${N_cores}}; 
export NUM_PROC

# RSYNCV="" the flag for rsync 
RSYNCV=""
if [ $V -ge 1 ]; then
    export RSYNCV="v" 
fi

# 
: ${TPL_MAIN_BRANCH:="master"} 
export TPL_MAIN_BRANCH

# 
: ${TPL_PKG:="thirdPartyLibs"} 
export TPL_PKG

: ${PKG:="GEOS"} ;
export PKG

# TPL_UPDATE:="1" => update TPL repos
: ${TPL_UPDATE:="1"} 
export TPL_UPDATE

# GEOS_UPDATE:="1" => update GEOS repos
: ${GEOS_UPDATE:="1"} 
export GEOS_UPDATE

# CVX_HOSTCONFIGS :="1" => add local host_configs/CVX to those from origin repos
: ${CVX_HOSTCONFIGS:="1"}; 
export CVX_HOSTCONFIGS

# TPL_REBASE:="1" => git pull --rebase TPL 
: ${TPL_REBASE:="0"} 
export TPL_REBASE

if [ $TPL_REBASE -eq 1 ]; then
    export TPL_PULL_OPTIONS=" --rebase "
fi
## TPL_RESET   '--type { HEAD~n | hash }' to reset before building TPL targets'

: ${GEOS_MAIN_BRANCH:="develop"} 
export GEOS_MAIN_BRANCH

# GEOS_REBASE:="1" => git pull --rebase GEOS 
: ${GEOS_REBASE:="0"} 
export GEOS_REBASE

if [ $GEOS_REBASE -eq 1 ]; then
    export GEOS_PULL_OPTIONS=" --rebase "
fi
## GEOS_RESET   '--type { HEAD~n | hash }' to reset before building GEOS targets'

: ${CLONE:="0"} 
export CLONE

# BUILD_ONLY : 0 = CLONE=1, clone and exit; 1 = build and do not clone 
: ${BUILD:="2"} 
if [ $BUILD -eq 1 ]; then
    BUILD_ONLY=1;
fi 

: ${BUILD_ONLY:="2"} 
if [ $BUILD_ONLY -eq 1 ]; then
    BUILD=1;
fi 
export CLONE BUILD_ONLY BUILD

: ${INIT_ONLY:="0"} 

if [[ $BUILD -eq 1 || $INIT_ONLY -eq 1 ]]; then
    PUBLISH_INIT=1;
fi 
export PUBLISH_INIT


# Shared libs default
: ${GEOS_BUILD_SHARED_LIBS:="ON"};
export GEOS_BUILD_SHARED_LIBS

# LD_GEOSX_CORE : 0 = no libgeosx_core.so in LD_LIBRARY_PATH; 1 = True
: ${LD_GEOSX_CORE:="0"};
export LD_GEOSX_CORE

# GPU GEOS workaround; default : NO workaround
: ${GEOS_GPU_WORKAROUND:="0"}; 
export GEOS_GPU_WORKAROUND

# Additional variables for the TPL Cmake configure command lines -DVAR=val
: ${TPL_CMAKE_VARS:=""};
export TPL_CMAKE_VARS

# Forcing these for GPU builds  -Xcompiler -DGEOS_USE_DEVICE  -Xcompiler -DCALC_FEM_SHAPE_IN_KERNEL 
# Additional variables for the GEOS Cmake configure command lines -DVAR=val
: ${GEOS_CMAKE_VARS:=""};
GEOS_CMAKE_VARS+=" -DENABLE_YAPF=OFF "
export GEOS_CMAKE_VARS

# Additional options for the host compilation command line, e.g., -DVAR=val
: ${GEOS_HOST_FLAGS_CLI:=""};
export GEOS_HOST_FLAGS_CLI

# Additional options for the host compilation command line, e.g., -DVAR=val
: ${GEOS_CUDA_FLAGS_CLI:=""};
export GEOS_CUDA_FLAGS_CLI

# TPL_BUILD_ONLY : 0 = configure and build all TPL packages, 'list ' = only rebuild those TPL packages in 'list'
: ${TPL_BUILD_ONLY:="0"}
export TPL_BUILD_ONLY 

# GEOS_REBUILD : 1 rebuild GEOS without configure step, 0 = configure and build all GEOS targets
: ${GEOS_REBUILD:="0"};
export GEOS_REBUILD

# GEOS_BUILD_ONLY : 0 = build all; 1 = only build geosx ; list of build targets 
: ${GEOS_BUILD_ONLY:="0"}; 
export GEOS_BUILD_ONLY

# GEOS_INSTALL : 0 = build and not install; 1 = build and install at install location ${iroot}/GEOS/XXX 
if [[ "$CHAP" == "1" ]]; then
    GEOS_INSTALL=1
fi

: ${GEOS_INSTALL:="1"}; 
export GEOS_INSTALL

: ${LCACHE:="1"}
export LCACHE 

: ${CTEST:="0"}
export CTEST 

if [ "${BUILD_ONLY}" == "0" ]; then
    export CTEST=0
fi 


# Python tools
# https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/index.html
: ${XML_TOOLS_BUILD:="0"}
export XML_TOOLS_BUILD

: ${GEOS_PYTHON_TOOLS:="${XML_TOOLS_BUILD}"}
export GEOS_PYTHON_TOOLS 
# make geosx_python_tools_clean
# make geosx_python_tools

: ${PYGEOS_BUILD:="0"}
export PYGEOS_BUILD

: ${GEOS_INTEGRATED_TESTS:="0"}
export GEOS_INTEGRATED_TESTS
# make ats_environment

## Release, Debug and RelWithDebInfo
: ${BUILD_TYPE:="RelWithDebInfo"} 
export BUILD_TYPE 

build_type=$(echo $BUILD_TYPE | tr '[:upper:]' '[:lower:]') 
export build_type

: ${ENABLE_TRILINOS:="OFF"}; 
export ENABLE_TRILINOS

: ${ENABLE_CALIPER_HYPRE:="OFF"}; 
export ENABLE_CALIPER_HYPRE

: ${ENABLE_CUDA_NVTOOLSEXT:="OFF"}; 
export ENABLE_CUDA_NVTOOLSEXT

# NOT used
: ${GEOSX_CONDA_ENV:=""}; 
export GEOSX_CONDA_ENV

: ${GEOSX_VENV:="geos-venv-${MPI}"}; 
export GEOSX_VENV

# Determine CUDA version in effect
export CUDA_VER=$([[ -x $(which nvcc 2>/dev/null) ]] && nvcc --version | gawk '$1=="Cuda" && $2=="compilation" {print $6}')


# Setting SIMPLE_PATHS to 0 to permanently disable this troublesome feature
export SIMPLE_PATHS="0"

if [ "${USER}" != "mtml" ]; then 
    # : ${SIMPLE_PATHS:="1"}; 
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
    # -------------------------------------------------------------------------- #
    : ${sroot:="/data/saet/${USER}/src"} ; 
    # : ${croot:="/mnt/resource/${USER}/src"} ; 
    : ${iroot:="/data/saet/hpcrnd/software"} ;
fi
: ${croot:="/dev/shm/${USER}/src"} ; 

## Bulding different branches (common or explicitly for TPL and/or GEOS)
# if [ -n "$BRANCH" ]; then 
#     : ${TPL_BRANCH:=${BRANCH}}
#     : ${GEOS_BRANCH:=${BRANCH}}
# fi 
if [ -n "${GEOS_BRANCH}" ]; then
    geos_branch=--$(echo $GEOS_BRANCH | sed  -e 's/ /_/g' -e 's|/|__|g')
    # geos_branch=--$(echo $GEOS_BRANCH | tr ' ' '_' | tr '/' '__')
    branch="${geos_branch}"
fi
if [ -n "${TPL_BRANCH}" ]; then
    tpl_branch=--$(echo $TPL_BRANCH | sed  -e 's/ /_/g' -e 's|/|__|g')
    # tpl_branch=--$(echo $TPL_BRANCH | tr ' ' '_' | tr '/' '__')
    branch+="${tpl_branch}"
fi

export GEOS_EFFECTIVE_BRANCH=${GEOS_BRANCH:=${GEOS_MAIN_BRANCH}}
export TPL_EFFECTIVE_BRANCH=${TPL_BRANCH:=${TPL_MAIN_BRANCH}}
export BRANCH TPL_BRANCH GEOS_BRANCH GEOS_EFFECTIVE_BRANCH

### ======================================================================= ####
## Last GEOS_BRANCH that builds
 ## Thu Sep 26 06:44:40 2024 +0200 : c225416b8ca7efcc08e13503f96dbc08e3c925f3
 ## Tue Oct 1 18:35:46 2024 +0200 : 0db85bed495f1b37bd1b64a4c0eef74d62717cc3
 ## Wed Oct 2 10:26:06 2024 -0500 : 678890e1dd8588e71667dc4afe1fcac0779ea4c1

### Does not build
 ## Thu Oct 3 05:23:41 2024 -0700 : 915d82f6b76cd1445a1ab115d919b879356a4671
 ## Mon Oct 7 15:13:22 2024 +0200 : 769dd04a5ceda609868e4675e84241278b1ea1eb
 ## Tue Nov 19 11:53:00 2024 -0600 : 6d285d2a8ecebb86bcc1459a3ccee179e4ce85ea
### ======================================================================= ####


# 
: ${projdir:=""} 

if [ "${projdir}" != "" ] ; then
    export projdir="-${projdir}"
fi

# projdir+=$branch

# -------------------------------------------------------------------------- #
# HOST_CONFIG_DIR is the file path prefix to host-configs directory 
# Clean up non-standard CONFIG_NAME and reconstruct "proper" full file path
export HOST_CONFIG_BASE=$(basename ${CONFIG_NAME} .cmake)
export HOST_CONFIG_DIR=$(readlink -f ${HOST_CONFIG_DIR})
export HOST_CONFIG_FILENAME="${HOST_CONFIG_DIR}/${HOST_CONFIG_BASE}.cmake"


# set umask 
umask 0022 

# -------------------------------------------------------------------------- #
# Create sroot or set it to pwd and manage source directories
# echo "## $0 PRE Note : sroot=$sroot "  ;
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

# Source path : under which TPL and GEOS are cloned
export srcpath="${sroot}/${PKG}"
if [ ! -d $srcpath ]; then		
    mkdir -p $srcpath ; RC=$?
    if [ $RC -ne 0 ]; then
	echo "## $0 Error : I cannot create $srcpath ! Exiting !  "  ;
	exit 20;
    fi
    echo "## $0 Note : Created srcpath=$srcpath "  ;
fi
export srcpath
echo "## $0 Note : Using srcpath=$srcpath  "  ;

# Where sources are actually cloned 
# export TPL_SOURCES_DIR="${srcpath}/${TPL_PKG}"
# export GEOSX_SOURCES_DIR="${srcpath}/${PKG}"
export TPL_SOURCES_DIR="${sroot}/${PKG}/${TPL_PKG}"
export GEOSX_SOURCES_DIR="${sroot}/${PKG}/${PKG}"

# Hack to retrieve latest version 
export geos_ver=$(gawk '{print $3}' ${GEOSX_SOURCES_DIR}/src/VERSION | sed 's/v//')
if [ -n "$geos_ver" ]; then
    export GEOS_VER=$geos_ver
fi
echo "## GEOS : GEOS_VER=$GEOS_VER"

# -------------------------------------------------------------------------- #
# Default environment path component names
: ${ARCH:="x86_64"} ;
: ${OS:="RHEL8"} ;
: ${PKG:="GEOS"} ;
: ${VER:="$GEOS_VER"} ;
# : ${VER:="0.2.0${projdir}"} ;
  STACK="${HOST_CONFIG_BASE}-${build_type}"
  BUILDDIR="build-${STACK}"
  INSTALLDIR="install-${STACK}"
  export GEOSX_VENV="${BUILDDIR}--$GEOSX_VENV" BUILDDIR INSTALLDIR
# -------------------------------------------------------------------------- #

# Manage and create local build cache directory
if [ ${LCACHE} -gt 0 ] ; then 
    [ ! -d $croot ] && mkdir -p $croot
    if [ ! -d $croot ]; then 
	echo "## $0 Warning : I cannot create local cache at $croot for builds; setting LCACHE=0 !" ;
	export LCACHE=0
    else
	free_cache=$(df -k $croot | gawk 'NR==2 {printf "%-l20f", $4/2^20}')
	free_cache_GiB=$(echo "$free_cache/1.0" | bc )
	if [ $free_cache_GiB -lt 30 ]; then
	    LCACHE=2
	    echo "## $0 Note : free_cache=$free_cache_GiB < 30 GiB; setting LCACHE=2 to clear $croot before builds " ;
	fi
	echo "## $0 Note : Using $croot to cache builds; free cache : $free_cache_GiB GiB ($free_cache GiB) " ;
    fi
fi

# Build directories
if [[ "${LCACHE}" == "0" ]] ; then
    broot=$sroot; 
else
    broot=$croot; 
fi
export TPL_BUILD_DIR=${broot}/${PKG}/${TPL_PKG}/$BUILDDIR
export GEOSX_BUILD_DIR=${broot}/${PKG}/${PKG}/$BUILDDIR

# Manage installation directories
if [[ "${CHAP}" == "1" ]]; then

    echo "## $0 Note : /chap public installation "  ;
    export iroot="/chap"
    echo "## $0 Note : Using iroot=$iroot "  ;

    isuffix_geosx="geos/${VER}${projdir}${geos_branch}${tpl_branch}/${ARCH}/${OS}/${INSTALLDIR}"
    prefix_geosx="${iroot}/${isuffix_geosx}" ; # [ ! -d $prefix_geosx ] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
    isuffix_tpl="TPL"
    prefix_tpl="${prefix_geosx}/${isuffix_tpl}" ; # [ ! -d $prefix_tpl ] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);
    mrootalli="${iroot}/geos/modulefiles"; # [ ! -d $mrootalli ] && mkdir -p $mrootalli; # mrootalli=$(readlink -f $mrootalli);
    initroot="${prefix_geosx}/.init"; # [ ! -d $initroot ] && mkdir -p $initroot; # initroot=$(readlink -f $initroot);
    initrooti="${iroot}/geos/.init"; # [ ! -d $initrooti ] && mkdir -p $initrooti; #  initrooti=$(readlink -f $initrooti);
    # initrootu="${HOME}/.init/${PKG}"; # [ ! -d $initrootu ] && mkdir -p $initrootu; # initrootu=$(readlink -f $initrootu);

else

    echo "## $0 Note : /data/saet/software installation "  ;
    if [[ ! -d $iroot ]]; then
	mkdir -p $iroot ; RC=$?
	if [ $RC -ne 0 ]; then
	    echo "## $0 Warning : I cannot create $iroot ! Setting iroot = $(pwd) "  ;
	    export iroot=$(pwd) 
	fi
    fi
    # iroot=$(readlink -f $iroot) ;
    echo "## $0 Note : Using iroot=$iroot "  ;

    isuffix_tpl="${ARCH}/${OS}/${PKG}TPL/${VER}${projdir}${tpl_branch}/${INSTALLDIR}"
    prefix_tpl="${iroot}/${isuffix_tpl}" ; # [ ! -d $prefix_tpl ] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);

    isuffix_geosx="${ARCH}/${OS}/${PKG}/${VER}${projdir}${geos_branch}${tpl_branch}/${INSTALLDIR}"
    # prefix_geosx="${iroot}/${suffix_geosx}" ; [ ! -d $prefix_geosx ] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
    
    if [[ "${GEOS}" == "1" && "${GEOS_INSTALL}" != "0" ]]; then 
	prefix_geosx="${iroot}/${isuffix_geosx}" ; 
	# [ ! -d $prefix_geosx ] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
    else
	prefix_geosx=$GEOSX_BUILD_DIR ; 
    fi
    mrootalli="${iroot}/modulefiles/${PKG}"; # [ ! -d $mrootalli ] && mkdir -p $mrootalli; # mrootalli=$(readlink -f $mrootalli);
    initroot="${prefix_geosx}/.init"; # [ ! -d $initroot ] && mkdir -p $initroot; # initroot=$(readlink -f $initroot);
    initrooti="${iroot}/.init/${PKG}"; # [ ! -d $initrooti ] && mkdir -p $initrooti; # initrooti=$(readlink -f $initrooti);
    # initrootu="${HOME}/.init/${PKG}"; # [ ! -d $initrootu ] && mkdir -p $initrootu; # initrootu=$(readlink -f $initrootu);
fi

# Installation directories
# Check and create TPL installation path
if [[ "${TPL}" == "1" && "${BUILD}" == "1" ]]; then
   [[ ! -d $prefix_tpl ]] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);
fi

# Check and create GEOS installation path
if [[ "${GEOS}" == "1" && "${BUILD}" == "1" && "${GEOS_INSTALL}" != "0" ]]; then 
   [[ ! -d $prefix_geosx ]] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
fi

# Allow to use existing TPL installation possibly from a common location by setting it to GEOSX_TPL_DIR
: ${GEOSX_TPL_DIR:=${prefix_tpl}} 
export GEOSX_TPL_DIR
export GEOS_TPL_DIR=$GEOSX_TPL_DIR

# GEOSX_DIR is computed
export GEOSX_DIR=${prefix_geosx}


if [[ ${V} -ge 1 ]]; then 
    echo "## $0 Note :  
#    prefix_geosx=$prefix_geosx ; isuffix_geosx=$isuffix_geosx 
#    prefix_tpl=$prefix_tpl ; isuffix_tpl=$isuffix_tpl"
fi

# -------------------------------------------------------------------------- #
modfile="${ARCH}-${OS}-${PKG}-${VER}${projdir}${geos_branch}${tpl_branch}-${STACK}"; 
# mrootallu="${HOME}/modulefiles/${PKG}"; # [ ! -d $mrootallu ] && mkdir -p $mrootallu; # mrootallu=$(readlink -f $mrootallu);

mroot="${prefix_geosx}/modulefiles"; # [ ! -d $mroot ] && mkdir -p $mroot;  # mroot=$(readlink -f $mroot);

initrcf="Init-${modfile}.rc";

# -------------------------------------------------------------------------- #
export broot croot sroot iroot mroot mrootalli initrooti isuffix_tpl isuffix_geosx prefix_geosx ARCH OS PKG VER STACK
# export broot croot sroot iroot mroot mrootallu mrootalli initrootu initrooti isuffix_tpl isuffix_geosx prefix_geosx ARCH OS PKG VER STACK


echo "# ==========================================================================================="
echo "# Chevron Automated TPL and GEOS clone-configure-build-install framework"
echo "#"
echo "#    $0 $@ "
echo "#"
echo "# Chevron Technology Center, Innovation and HPC R&D "
echo "#"
echo "# --------------- Configuration Status -------------------------------------------------------"
echo "#"
echo "# CONFIG_NAME          = $CONFIG_NAME "
echo "#"
echo "# HOST_CONFIG_BASE     = $HOST_CONFIG_BASE"
echo "# HOST_CONFIG_DIR      = $HOST_CONFIG_DIR"
echo "# HOST_CONFIG_FILENAME = $HOST_CONFIG_FILENAME"
echo "# TPL_SOURCES_DIR      = ${TPL_SOURCES_DIR} ; where TPL sources are cloned"
echo "# GEOSX_SOURCES_DIR    = ${GEOSX_SOURCES_DIR} ; where GEOS sources are cloned"
echo "# TPL_BUILD_DIR        = ${TPL_BUILD_DIR} ; where TPL will be built"
echo "# GEOSX_BUILD_DIR      = ${GEOSX_BUILD_DIR} ; where GEOS will be built"
echo "# GEOSX_TPL_DIR        = ${GEOSX_TPL_DIR} ; if set, it will be used to let GEOS targets link against TPL at this existing library location "
echo "# GEOS_TPL_DIR         = ${GEOS_TPL_DIR} ; same as GEOSX_TPL_DIR "
echo "# GEOSX_DIR            = ${GEOSX_DIR} ; where GEOS will be installed"
echo "# GEOSX_VENV           = ${GEOSX_VENV} ; name of virtual environment for Python packages"
echo "# GEOSX_INSTALL_VENV_DIR it will appear later before the virtual environment for the Python tools gets built "
echo "# CVX_PYTHON3          = ${CVX_PYTHON3} "

if [ ${V} -gt 0 ]; then

    source ${scrdir}/GEOS-quick-help.sh 

    echo "# --------------------------------------------------------------------------------------------"
    env | egrep -i 'slurm(.*)=' | sort | gawk '{print "# ", $0 }'
    echo "# --------------------------------------------------------------------------------------------"
    echo "# sroot           = ${sroot}" 
    echo "# croot           = ${croot}" 
    echo "# broot           = ${broot}" 
    echo "# iroot           = ${iroot}" 
    echo "# prefix_tpl      = ${prefix_tpl} "
    echo "# prefix_geosx    = ${prefix_geosx} "
    echo "# projdir         = ${projdir}" 
    echo "# modfile         = ${modfile}" 
    echo "# srcpath         = ${srcpath} ; the file directory root where GEOS and TPL will be cloned" 
    echo "# scrpath dir     = ${scrpath} ; the file directory root these scripts reside on " 
    echo "# Script dir/name = ${scrdir}/$(basename $0) ; full pathname of this script"
    echo "# BUILDDIR         = ${BUILDDIR}" 
    echo "# INSTALLDIR      = ${INSTALLDIR}" 
    echo "# LCACHE          = ${LCACHE} : 0 = no cache; 1 = use cache; 2 = clear cache before building; 3 = copy GEOS cached build to \$sroot; 4 = copy entire cache back" 
    echo "# CTEST           = ${CTEST} : 0 = don't run tests; 1 = run tests " 
    if [ ${V} -gt 1 ]; then 
	echo "# --------------------------------------------------------------------------------------------"
	env | egrep 'comp|MODULE|mpi|GEOS|CUDA|(ARCH|OS|PKG|VER|STACK|(I_|O)MPI_.*)=' | sort | gawk '{print "# ", $0 }'
	printf "# mpicc           = $(which mpicc)  \t[ Intel : $(which mpiicc 2>/dev/null) ] \n"
	printf "# mpicxx          = $(which mpicxx) \t[ Intel : $(which mpiicpc 2>/dev/null) ] \n"
	printf "# mpifort         = $(which mpifort)\t[ Intel : $(which mpiifort 2>/dev/null) ] \n"
    fi
fi
echo "# --------------------------------------------------------------------------------------------"

# 
if [ $# -lt 1 ]; then 
    echo "## $0 Error : Backend script may only be sourced by the configuration scripts ! Exiting ... "  ;
    exit 40; 
fi 

# At this stage $srcpath is the location to receive the thriParyLibs and GEOSX repos
export DT_GEOSX=$(date +%F-%H%M%S)

TS_0=$(date "+ %s.%N")
TS_TPL_0="0.0"
TS_TPL_1="0.0"
TS_CTEST_0="0.0"
TS_CTEST_1="0.0"
    
TS_GEOS_0="0.0"
TS_GEOS_1="0.0"

echo "## ${HOST_CONFIG_BASE} Start : ${DT_GEOSX} ${TS_0} "
# cd to root where TPL and GEOS are cloned
cd $srcpath
echo "## Now at $(pwd) ( must == $srcpath )"

if [ "${INIT_ONLY}" != "1" ]; then 

    git-lfs install  || RCerror $0 $? "(GEOS) Cannot git install "

    TS_0=$(date "+ %s.%N")

    if [ "${GEOS_REBUILD}" != "1"   ] ; then 
	# Check and clone, update, configure, build and install TPL and GEOS
	if [ "${CLONE}" == "1" ]; then
	    echo "## GEOS : Cloning $gitGEOS --> $(pwd)/$PKG "
	    git clone ${gitGEOS} $PKG ; RC=$?
 	    if [[ $RC != 0 ]]; then
		echo "## GEOS Error RC=$RC : Cloning $gitGEOS --> $(pwd)/$PKG failed! Exiting ,,,"
		exit $RC;
	    fi
	fi

	if [[ "${GEOS_UPDATE}" == "1" && "${CLONE}" != "1" ]] ; then 
	    cd ${PKG}
	    echo "## GEOS : Updating ${PKG} at $(pwd) (also update the GEOS/host-configs before configuring and building TPL) "
	    echo "## Origin : ${gitGEOS}"

	    git checkout $GEOS_MAIN_BRANCH || RCerror $0 $? "(GEOS) Cannot checkout $GEOS_MAIN_BRANCH"

	    git pull origin $GEOS_PULL_OPTIONS || RCerror $0 $? "(GEOS) Cannot git pull origin $GEOS_PULL_OPTIONS "

	    git lfs pull origin || RCerror $0 $? "(GEOS) Cannot git lfs pull origin "

	    git submodule init

	    # Building integratedTests now mandatory
	    # if [ "$GEOS_INTEGRATED_TESTS" == "0" ] ; then 
		# git submodule deinit integratedTests
	    # fi

	    git submodule update

	    ## Detailed submodule updates 
	    git submodule update --init src/cmake/blt
	    git submodule update --init src/coreComponents/LvArray
	    git submodule update --init src/coreComponents/fileIO/coupling/hdf5_interface
	    git submodule update --init src/coreComponents/constitutive/PVTPackage
	    # git submodule update --init scripts/uberenv
	    
	    # Populate $srcpath/$PKG/host-configs/CVX with CVX specific host_configs 
	    if [[ $CVX_HOSTCONFIGS == 1 ]]; then 
		echo "## GEOS : Populating CVX specific $scrdir/$PKG/host-configs/CVX  --> $srcpath/$PKG/host-configs/  "
		rsync -a${RSYNCV} $scrdir/$PKG/host-configs/CVX $srcpath/$PKG/host-configs/
		# echo "## GEOS : Populating CVX specific $scrdir/$PKG/host-configs/CVX  --> $srcpath/$PKG/host-configs/  "
		# rsync -a${RSYNCV} $scrdir/$PKG/host-configs/CVX $srcpath/$PKG/host-configs/
		ls -l $srcpath/$PKG/host-configs/tpls.cmake $srcpath/$PKG/host-configs/CVX/tpls.cmake | gawk '{print "#   ", $0}'
		echo "## GEOS : Updating local host-configs with latest public GihHub : $srcpath/$PKG/host-configs --> $scrdir/$PKG/ excluding our's at ./CVX/* "
		rsync -a${RSYNCV} $srcpath/$PKG/host-configs --exclude '*CVX/*' $scrdir/$PKG/
	    else
		echo "## GEOS : Did NOT populate CVX specific $scrdir/$PKG/host-configs/CVX  --> $srcpath/$PKG/host-configs/  "
	    fi
	fi
    else
	echo "## GEOS Note : NOT updating ${PKG} at $(pwd) ( GEOS_REBUILD = ${GEOS_REBUILD}) "
    fi
    
    ## TPL clone, configure and build-install tasks -------------------------------------------------------------------------
    if [ "${TPL}" == "1" ] ; then 
	TS_TPL_0=$(date "+ %s.%N")
	# cd to root where TPL and GEOS are cloned
	cd $srcpath
	echo "## TPL : at $(pwd) ( must == $srcpath )"

	if [ "${CLONE}" == "1" ]; then
	    echo "## TPL : Cloning ${gitTPL} --> ${TPL_PKG}"
	    git clone ${gitTPL} ${TPL_PKG} ; RC=$?
 	    if [[ $RC != 0 ]]; then
		echo "## TPL Error : Cloning ${gitTPL} --> ${TPL_PKG} failed! Exiting ,,,"
		exit $RC;
	    fi
	    # https://github.com/GEOS/thirdPartyLibs.git
	fi

	if [[ "${TPL_UPDATE}" == "1" && "${CLONE}" != "1" ]] ; then 
	    # Update the TPL packages 
	    cd ${TPL_SOURCES_DIR}
	    echo "## TPL : Updating ${TPL_PKG} at $(pwd) ( should== $TPL_SOURCES_DIR ) "
	    echo "## TPL Origin : ${gitTPL}"

	    git checkout $TPL_MAIN_BRANCH || RCerror $0 $? "(TPL) Cannot git checkout $TPL_MAIN_BRANCH"

	    git-lfs install || RCerror $0 $? "(TPL) Cannot git-lfs install"

	    git pull origin $TPL_PULL_OPTIONS || RCerror $0 $? "(TPL) Cannot git pull origin $TPL_PULL_OPTIONS "

	    git-lfs pull origin || RCerror $0 $? "(TPL) Cannot git-lfs pull origin "

	    git submodule init
	    git submodule update
	fi

	## Checkout required TPL branch to work with
	if [ -n "${TPL_BRANCH}" ]; then
	    echo "## TPL : Checking out branch=${TPL_BRANCH} "
	    cd ${TPL_SOURCES_DIR}
	    echo "## TPL : at $(pwd) ( should== $TPL_SOURCES_DIR ) "
	    git checkout $TPL_BRANCH ; RC=$?
	    if [ $RC -ne 0 ]; then
		echo "## $0 TPL ERROR : There is no TPL branch $TPL_BRANCH at $(pwd) ! Exititng ... "
		exit $RC;
	    fi

	    if [[ "${TPL_UPDATE}" == "1" && "${CLONE}" != "1" ]] ; then 
		echo "## TPL : Updating branch (git pull branch=${TPL_BRANCH}) "
		git pull || RCerror $0 $? "(TPL) Cannot git pull "
		git-lfs pull || RCerror $0 $? "(TPL) Cannot git-lfs pull "
	    
		git submodule init
		git submodule update
	    fi
	fi

	if [ "${BUILD_ONLY}" == "1" ]; then
	    ## Configure and Build 

	    if [ ${LCACHE} -gt 0 ]; then
		if [ ${LCACHE} -eq 2 ]; then
		    echo "## TPL : Clearing $croot/$PKG/${TPL_PKG} "
		    \rm -fr $croot/$PKG/${TPL_PKG}
		fi

		echo "## TPL : Syncing $sroot/$PKG/${TPL_PKG} --> $croot/${PKG} "
		rsync --exclude 'thirdPartyLibs/*build*' --exclude 'thirdPartyLibs/*install*' -a${RSYNCV} $sroot/$PKG/${TPL_PKG} $croot/$PKG

		echo "## TPL : Syncing $scrdir/$PKG/host-configs --> $croot/$PKG/$PKG "
		rsync -a${RSYNCV} $scrdir/$PKG/host-configs $croot/$PKG/$PKG

		cd $croot/$PKG/${TPL_PKG} 
		echo "## Now at $(pwd) "
	    fi

	    # At either TPL_SOURCES_DIR or $croot/$PKG/${TPL_PKG}
	    # BUILDDIR="build-${CONFIG_NAME}-${build_type}"
	    if [ ! -d ${BUILDDIR} ]; then
		echo "## TPL : Creating ${BUILDDIR} at $(pwd) "
		mkdir -p ${BUILDDIR}; 
	    fi

	    # Which parts of TPL to build ; 0 ; all; or, target list
	    if [ "${TPL_BUILD_ONLY}" == "0" ]; then 

		if [ -n "${TPL_RESET}" ]; then 
		    echo "## TPL : git reset ${TPL_RESET} "
		     git reset ${TPL_RESET} 
		fi

		echo "## TPL : Configuring thirdPartyLibs at $(pwd) (with SCOTCH_NUM_PROC=1)"
		# includes workaround -DSCOTCH_NUM_PROC=1
		$ECHO $CVX_PYTHON3 scripts/config-build.py "$TPL_CMAKE_VARS" -hc $HOST_CONFIG_FILENAME -DSCOTCH_NUM_PROC=1 \
		    -bt ${BUILD_TYPE} -DGEOSX_TPL_DIR=$GEOSX_TPL_DIR -ip $GEOSX_TPL_DIR -DNUM_PROC=${NUM_PROC} |& tee Config-TPL-${DT_GEOSX}.log  
 
		cd ./${BUILDDIR}
		export TPL_BUILDDIR_PATH=$(readlink -f $(pwd))
		echo "## TPL : Building thirdPartyLibs at TPL_BUILD_DIR=$TPL_BUILD_DIR ( TPL_BUILDDIR_PATH=$TPL_BUILDDIR_PATH == $TPL_BUILD_DIR ) "
		make -k |& tee Make-${DT_GEOSX}.log 
	    else 
		# Only rebuild list of TPL targets
		cd ./${BUILDDIR}
		export TPL_BUILDDIR_PATH=$(readlink -f $(pwd))
		echo "## TPL : Re-Building thirdPartyLibs at TPL_BUILD_DIR=$TPL_BUILD_DIR ( TPL_BUILDDIR_PATH=$TPL_BUILDDIR_PATH == $TPL_BUILD_DIR ) "
		VERBOSE=1 make ${TPL_BUILD_ONLY} -k |& tee Make-TPL_BUILD_ONLY-${DT_GEOSX}.log 
	    fi
	    ln -s $GEOSX_TPL_DIR

	 build_configuration="TPL-build-configuration--$(hostname)-${DT_GEOSX}.txt"
	 echo "## TPL built on $(hostname) in $(pwd) at ${DT_GEOSX} by ${USER} 
-------------------------------------------------------------------------------------------------------
$(cat cmake_cmd) 
-------------------------------------------------------------------------------------------------------
# $(module list)
-------------------------------------------------------------------------------------------------------
#  sroot=${sroot} 
#  iroot=${iroot} 
#  croot=${croot}  
$(env | egrep 'root|fix|Python3|PYTH|ENV|ENABLE(.*)|TPL|GEOS|CONF|CACHE|CUDA|ARCH|OS|PKG|VER|STACK|HPCX|MPI|STACK' | sort | gawk '{print "# ", $0 }')" |& tee -a $GEOSX_TPL_DIR/${build_configuration}
	fi

	TS_TPL_1=$(date "+ %s.%N")
    fi 

    ## GEOS configure, build and install tasks  ---------------------------------------------------------------------------------
    cd $srcpath 
    echo "## Now at $(pwd) ( should == $sroot/${PKG})"

    if [ "${GEOS}" == "1" ] ; then 
	TS_GEOS_0=$(date "+ %s.%N")

	if [ "${BUILD_ONLY}" == "1" ]; then
	    ## Configure and Build 
	    echo "## GEOS : Preparing to build ${PKG}"
	    cd $GEOSX_SOURCES_DIR
	    
	    ## Checkout required GEOS branch to work with
	    if [ -n "${GEOS_BRANCH}" ]; then

		echo "## GEOS : Checking out branch=${GEOS_BRANCH} at $(pwd) "
		git checkout $GEOS_BRANCH ; RC=$?
		if [ $RC -ne 0 ]; then
		    echo "## $0 GEOS ERROR : There is no GEOS branch $GEOS_BRANCH at $(pwd) ! Exititng ... "
		    exit $RC;
		fi
		if [ "${GEOS_UPDATE}" == "1" ] ; then 
		    echo "## GEOS : updating branch=${GEOS_BRANCH} (git pull) at $(pwd) "
		    git pull || RCerror $0 $? "(GEOS) Cannot git pull"
		    git-lfs pull || RCerror $0 $? "(GEOS) Cannot git-lfs pull"
		    git submodule init || RCerror $0 $? "(GEOS) Cannot submodule init "
		    git submodule update || RCerror $0 $? "(GEOS) Cannot submodule update "

		else
		    echo "## GEOS : NOT updating branch=${GEOS_BRANCH} from origin at $(pwd) "
		fi
	    fi

	    if [ -n "${GEOS_RESET}" ]; then 
		echo "## GEOS : git reset ${GEOS_RESET} at $(pwd) "
		git reset ${GEOS_RESET} 
	    fi

	    if [ ${LCACHE} -eq 2 ]; then
		echo "## GEOS : Clearing $croot/$PKG/$PKG "
		\rm -fr $croot/$PKG/$PKG
	    fi

	    if [ ${LCACHE} -gt 0 ]; then
		echo "## GEOS : Syncing $sroot/$PKG/$PKG --> $croot/$PKG "
		if [[ "$GEOS_INTEGRATED_TESTS" == "1"  ]]; then 
		    ## rsync  --exclude "$PKG/src/docs/*" --exclude "$PKG/*build*" --exclude "$PKG/*install*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		    rsync --exclude "$PKG/*build*" --exclude "$PKG/*install*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		else
		    ## rsync  --exclude "$PKG/src/docs/*" --exclude "$PKG/*build*" --exclude "$PKG/*install*" --exclude "$PKG/*inputFiles*"  --exclude "$PKG/*integratedTests*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		    rsync --exclude "$PKG/*build*" --exclude "$PKG/*install*" --exclude "$PKG/*inputFiles*"  --exclude "$PKG/*integratedTests*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		fi
		cd $croot/$PKG/$PKG 
	    fi

	    ## Now $sroot/$PKG/$PKG or at $croot/$PKG/$PKG
	    ## GEOS Configuration step
	    echo "## GEOS : Now at $(pwd) "

	    if [ ! -d ${BUILDDIR} ]; then
		echo "## GEOS : Creating BUILDDIR=${BUILDDIR} at $(pwd) "
		mkdir -p ${BUILDDIR};
	    fi
	    echo "## GEOS : BUILDDIR=${BUILDDIR} at $(pwd) "
	    ## MikeT: Check if the GEOSX_BUILD_DIR is the same as the one already set
	    echo "## GEOS : GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (Pre)"
	    export GEOSX_BUILD_DIR=$(readlink -f ${BUILDDIR})
 	    echo "## GEOS : GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (Post)"
	    
	    # GEOS_REBUILD == 1 does not install/rebuild ATS or Python Tools
	    if [ "${GEOS_REBUILD}" != "1"   ] ; then 

		if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o "${PYGEOS_BUILD}" != "0" ]; then
		    # Prepare virtual environment for Python tools and Designate compiler for mpi4py 
		    # GEOSX_BUILD_VENV

		    # Setup virtual environment as ../$BUILDDIR--$GEOSX_VENV a peer dir to $BUILDDIR
		    export GEOSX_BUILD_VENV_DIR=${GEOSX_TPL_DIR}/$GEOSX_VENV
		    export GEOSX_INSTALL_VENV_DIR=${GEOSX_BUILD_VENV_DIR}
		    
		    echo "## GEOS : Preparing GEOS Python virtual environment $GEOSX_VENV  (at $(pwd)) [GEOSX_BUILD_VENV_DIR=$GEOSX_BUILD_VENV_DIR]"
		    echo "##        Installing GEOS Python virtual environment at $GEOSX_BUILD_VENV_DIR [GEOSX_INSTALL_VENV_DIR=$GEOSX_INSTALL_VENV_DIR]"
		    # $CVX_PYTHON3 -m pip install venv
		    $CVX_PYTHON3 -m venv $GEOSX_BUILD_VENV_DIR
		    if [ -r $GEOSX_BUILD_VENV_DIR/bin/activate ]; then 
			source $GEOSX_BUILD_VENV_DIR/bin/activate
			RC=$?
			$CVX_PYTHON3 -m pip install --upgrade pip
			if [ $RC -eq 0 ]; then 
			    echo "## GEOS : Activated GEOS Python virtual environment $GEOSX_VENV at $GEOSX_BUILD_VENV_DIR"
			    echo "##        Using Python = $CVX_PYTHON3 "
			    export CVX_PYTHON3=$(which python)
			    export Python3_ROOT_DIR=$GEOSX_BUILD_VENV_DIR
			    export Python3_EXECUTABLE=$(basename ${CVX_PYTHON3})
			    export GEOS_CMAKE_VARS=" -DPython3_ROOT_DIR=$Python3_ROOT_DIR -DPython3_EXECUTABLE=$Python3_EXECUTABLE ${GEOS_CMAKE_VARS} "
			    if [ $V -gt 0 ]; then 
				echo "##  Python3_ROOT_DIR = $Python3_ROOT_DIR "
				echo "##  Python3_EXECUTABLE = $Python3_EXECUTABLE "
				echo "##  CVX_PYTHON3 = $CVX_PYTHON3 "
				echo "##  GEOS_CMAKE_VARS = \"${GEOS_CMAKE_VARS}\" "
			    fi
			else
			    echo "## GEOS ERROR : Failed activating GEOS Python virtual environment $GEOSX_VENV at $GEOSX_BUILD_VENV_DIR"
			    echo "## User Python packages will be installed at ${HOME}/.local/lib and binaries at ${HOME}/.local/bin"
  			fi

		    else
			echo "## GEOS Python ERROR : Could not activate GEOS Python virtual environment $GEOSX_VENV at $GEOSX_BUILD_VENV_DIR ! "
			echo "## User Python packages will be installed at ${HOME}/.local/lib and binaries at ${HOME}/.local/bin"
			# exit 10;
		    fi
		fi
		
		echo "## GEOS : Configuring GEOS at $(pwd) ( should== $GEOSX_BUILD_DIR ) using Python=$CVX_PYTHON3"
		
		$ECHO $CVX_PYTHON3 scripts/config-build.py "$GEOS_CMAKE_VARS" -hc $HOST_CONFIG_FILENAME -bt ${BUILD_TYPE} -ip $GEOSX_DIR -DGEOSX_DIR=$GEOSX_DIR -DGEOSX_TPL_DIR=$GEOSX_TPL_DIR |& tee Config-${DT_GEOSX}.log 

	    else
		echo "## GEOS Note : At this state a 'configure' step must have already taken place for ${CONFIG_NAME} "
	    fi
	    
	    if [ "${GEOS_BUILD_ONLY}" == "1" ]; then
		export GEOS_BUILD_ONLY="geosx"
	    elif [ "${GEOS_BUILD_ONLY}" == "0" ]; then
		export GEOS_BUILD_ONLY="all"
	    fi

	    ## GEOS Build step
	    echo "## GEOS : Building GEOS targets ${GEOS_BUILD_ONLY}  "
	    cd $GEOSX_BUILD_DIR
 	    echo "## GEOS : Now at $(pwd) ( should== $GEOSX_BUILD_DIR )"
	    make -j ${NUM_PROC} -k ${GEOS_BUILD_ONLY} |& tee Make-build-${DT_GEOSX}.log 

	    if [ "${GEOS_INTEGRATED_TESTS}" != "0" ]; then
		export GEOS_PYTHON_TOOLS=1
	    fi

	    # Prepare virtual environment for Python tools and Designate compiler for mpi4py 
	    # GEOSX_BUILD_VENV_DIR
	    
	    # Prepare virtual environment for Python tools and Designate compiler for mpi4py 
	    if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o  "${PYGEOS_BUILD}" != "0" ]; then
		# pip install ../src/coreComponents/python/* 
		export MPICC=$GEOSX_MPICC
		# echo "## GEOS : Removing existing python package mpi4py from $GEOSX_BUILD_VENV_DIR (at $(pwd) [Python : $CVX_PYTHON3])"
		# 1) $CVX_PYTHON3 -m pip cache remove mpi4py
		# 2) VERBOSE=1 env MPICC=$GEOSX_MPICC $CVX_PYTHON3 -m pip uninstall mpi4py
		echo "## GEOS : Building new python package mpi4py in $GEOSX_BUILD_VENV_DIR (at $(pwd) [Python : $CVX_PYTHON3])"
		echo "## GEOS : MPICC= $MPICC)"
		# 3) env LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64:/usr/lib" MPICC=$GEOSX_MPICC $CVX_PYTHON3 -m pip install --upgrade --no-cache-dir mpi4py
		env LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64:/usr/lib" MPICC=$GEOSX_MPICC $CVX_PYTHON3 -m pip install --upgrade --no-cache-dir mpi4py
		# export GEOSX_CONDA_ENV="${GEOSX_DIR}/Conda_Env"

		echo "## GEOS : Building GEOS Python/XML tools at $(pwd); All executables are installed in ${GEOSX_BUILD_VENV_DIR}/bin "
		export PATH="${GEOSX_BUILD_VENV_DIR}/bin:${PATH}"
		# export PATH="${HOME}/.local/bin:${PATH}"
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make geosx_python_tools_clean  |& tee -a Make-build-${DT_GEOSX}.log
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make geosx_python_tools |& tee -a Make-build-${DT_GEOSX}.log
	    fi

	    if [ "${GEOS_INTEGRATED_TESTS}" != "0" ]; then
		# export PATH="${HOME}/.local/bin:${PATH}"
		echo "## GEOS : Building intregratedTests (ats_environment) at $GEOSX_BUILD_VENV_DIR "
		# export LD_LIBRARY_PATH="/usr/lib64:/usr/lib:${LD_LIBRARY_PATH}"
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make ats_environment |& tee -a Make-build-${DT_GEOSX}.log
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make ats_clean |& tee -a Make-build-${DT_GEOSX}.log
		###  OMP_NUM_THREADS=1 ./geos_ats.sh --machine openmpi --ats openmpi_args=--oversubscribe --ats openmpi_procspernode=4 --ats openmpi_maxprocs=8 --ats openmpi_args=--report-bindings --ats openmpi_args="--bind-to none" --ats openmpi_mpirun $(which mpirun) --ats openmpi_install $HPCX_MPI_DIR
	    fi
 
	    ## Installation Phase
	    if [ "${GEOS_INSTALL}" == "1" ]; then 
		echo "## GEOS : Installing GEOS in $GEOSX_DIR "
		ln -s $GEOSX_DIR
		
		[[ -r  $GEOSX_TPL_DIR/${build_configuration} ]] && cp $GEOSX_TPL_DIR/${build_configuration} $GEOSX_DIR

		if [[ -x $GEOSX_DIR/bin/geosx ]]; then
		    echo "## GEOS : Removing old $GEOSX_DIR/bin/geosx "
		    \rm -f $GEOSX_DIR/bin/geosx
		fi
		if [[ -x $GEOSX_DIR/bin/geos ]]; then
		    echo "## GEOS : Removing old $GEOSX_DIR/bin/geos "
		    \rm -f $GEOSX_DIR/bin/geos
		fi
		     
		make -k install |& tee -a Make-install-${DT_GEOSX}.log 

		## Check for installation issues
		if [ ! -x $GEOSX_DIR/bin/geosx -o "${GEOS_GPU_WORKAROUND}" == "1" ]; then
		    echo "## GEOS Note : Install workaround "    |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./include $GEOSX_DIR  |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./bin $GEOSX_DIR      |& tee -a Make-install-${DT_GEOSX}.log 
		    cp  ./bin/geosx $GEOSX_DIR/bin/geosx    |& tee -a Make-install-${DT_GEOSX}.log 
		    cp  ./bin/geos $GEOSX_DIR/bin/geos    |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./tests/* $GEOSX_DIR/bin |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./lib $GEOSX_DIR      |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./share $GEOSX_DIR    |& tee -a Make-install-${DT_GEOSX}.log 
		    mkdir -p $GEOSX_DIR/lib64               |& tee -a Make-install-${DT_GEOSX}.log 
		    cp ./lib/libbenchmark* $GEOSX_DIR/lib64 |& tee -a Make-install-${DT_GEOSX}.log 
		    $GEOSX_DIR/bin/geosx -s $GEOSX_DIR/this-schema.xsd |& tee -a Make-install-${DT_GEOSX}.log 
		fi
		if [ ! -x $GEOSX_DIR/bin/geos ]; then
		    \cp -f $GEOSX_DIR/bin/geosx $GEOSX_DIR/bin/geos   |& tee -a Make-install-${DT_GEOSX}.log 
		fi
		# Propagate Python's virtual environment properly to install location
		# ATS may only be used at the build location !
		if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o  "${PYGEOS_BUILD}" != "0" ]; then
		    echo "## GEOS : Propagating Python Tools of GEOS to $GEOSX_DIR " |& tee -a Make-install-${DT_GEOSX}.log
		    echo "## GEOS : python virtual environment $GEOSX_VENV installed at GEOSX_INSTALL_VENV_DIR "  |& tee -a Make-install-${DT_GEOSX}.log
		    set -x
		    \rm -f $GEOSX_DIR/bin/geosx_preprocessed 
		    \cp -f ${GEOSX_SOURCES_DIR}/scripts/automatic_xml_preprocess.sh $GEOSX_DIR/bin/geosx_preprocessed  |& tee -a Make-install-${DT_GEOSX}.log
		    \rm -f $GEOSX_DIR/bin/pygeosx_preprocess.py 
		    \cp -f ${GEOSX_SOURCES_DIR}/scripts/pygeosx_preprocess.py  $GEOSX_DIR/bin/pygeosx_preprocess.py  |& tee -a Make-install-${DT_GEOSX}.log
		    chmod a+x $GEOSX_DIR/bin/pygeosx_preprocess.py
		    # \cp -f $(readlink -f ./bin/geosx_preprocessed) $GEOSX_DIR/bin/geosx_preprocessed  |& tee -a Make-install-${DT_GEOSX}.log
		    # \cp -f $(readlink -f ./bin/pygeosx_preprocess.py) $GEOSX_DIR/bin  |& tee -a Make-install-${DT_GEOSX}.log
		    set +x 
		    echo "## GEOS NOTE : ATS may only be used out of the build location in $GEOSX_BUILD_VENV_DIR !"
		    echo "##             Other Python tools could be used out of the install location in GEOSX_INSTALL_VENV_DIR or $GEOSX_DIR/bin !"
		fi
		## ATS is not relocatable and may only be usd at its installation location (peer dir where GEOS is built)
		# if [ "${GEOS_INTEGRATED_TESTS}" != "0" ]; then
		#     echo "## GEOS : Preparing intregratedTests (ats_environment) at $GEOSX_INSTALL_VENV_DIR "
		#     rsync -a${RSYNCV} integratedTests $GEOSX_DIR
		#     if [ -h $GEOSX_DIR/integratedTests/integratedTests ]; then
		# 	\rm  $GEOSX_DIR/integratedTests/integratedTests
		# 	echo "## GEOS : Linking $GEOSX_DIR/integratedTests/integratedTests -> $GEOSX_DIR/integratedTests/workingDir"
		# 	ln -fs $GEOSX_DIR/integratedTests/workingDir $GEOSX_DIR/integratedTests/integratedTests 
		#     fi
		# fi
	    else
		# NO installation of GEOS to GEOSX_DIR, adjust prefix_geosx to point to the build directory
		# if [[ "$LCACHE" == "1" ||  "$LCACHE" == "2" ]] ; then 
		#     export prefix_geosx=$(pwd)
		#     export GEOSX_DIR=${GEOSX_BUILD_DIR}
		# else
		#     export prefix_geosx=$(pwd)
		#     export GEOSX_DIR=${GEOSX_BUILD_DIR}
		# fi
		echo "## GEOS Note: NO installation, setting prefix_geosx=$prefix_geosx (the build directory: $(pwd))"
		echo "##            LCACHE = $LCACHE"
		echo "##            prefix_geosx = $prefix_geosx"
		echo "##            GEOSX_DIR = $GEOSX_DIR"
		echo "##            GEOSX_INSTALL_VENV_DIR = $GEOSX_INSTALL_VENV_DIR"
	    fi
	    build_configuration="GEOS-build-configuration--$(hostname)-${DT_GEOSX}.txt"
	    echo "## $PKG built on $(hostname) with developer stack $STACK in $(pwd) at ${DT_GEOSX} by ${USER} 
-------------------------------------------------------------------------------------------------------
# $0 $@
-------------------------------------------------------------------------------------------------------
$(cat cmake_cmd) 
-------------------------------------------------------------------------------------------------------
$(module list)
-------------------------------------------------------------------------------------------------------
#  sroot=${sroot} 
#  iroot=${iroot} 
#  croot=${croot} 
-------------------------------------------------------------------------------------------------------
$(env | egrep 'root|fix|Python3|PYTH|ENV|ENABLE(.*)|TPL|GEOS|CONF|CACHE|CUDA|ARCH|OS|PKG|VER|STACK|HPCX|MPI|STACK' | sort | gawk '{print "# ", $0 }') " |& tee -a $GEOSX_DIR/${build_configuration}
	fi

	TS_GEOS_1=$(date "+ %s.%N")

    fi 

    ## Run ctest -V  

    if [ ${CTEST} -gt 0 ]; then
	TS_CTEST_0=$(date "+ %s.%N")
	pushd $GEOSX_BUILD_DIR

	# if [ ${LCACHE} -gt 0 ]; then
	#     pushd $croot/$PKG/$PKG/$BUILDDIR
	# else
	#     pushd $sroot/$PKG/$PKG/$BUILDDIR
	# fi
	echo "## GEOS : Running ctest -V on $BUILDDIR at $GEOSX_BUILD_DIR. If there is no $BUILDDIR, this will fail."
	ctest -V |& tee CTEST-${CONFIG_NAME}_${DT_GEOSX}.log 

	if [ ${LCACHE} -gt 0 ]; then
	    echo "## GEOS : Syncing $croot/$PKG/$PKG/$BUILDDIR/Testing_${DT_GEOSX} --> $sroot/$PKG/$PKG/$BUILDDIR"
	    mkdir -p $sroot/$PKG/$PKG/$BUILDDIR
	    mv $croot/$PKG/$PKG/$BUILDDIR/Testing $croot/$PKG/$PKG/$BUILDDIR/Testing_${DT_GEOSX}
	    rsync -a${RSYNCV} $croot/$PKG/$PKG/$BUILDDIR/Testing_${DT_GEOSX} $sroot/$PKG/$PKG/$BUILDDIR    
	else
	    mv $sroot/$PKG/$PKG/$BUILDDIR/Testing $sroot/$PKG/$PKG/$BUILDDIR/Testing_${DT_GEOSX}
	fi 
	TS_CTEST_1=$(date "+ %s.%N")
    fi
    
    echo "## ${CONFIG_NAME} End  : $(date +%F-%H%M%S) $(date +%s.%N) "

    # 
    # set -x 
    TS_1=$(date "+ %s.%N")

    statsf=$GEOSX_DIR/${build_configuration} 
	
    dTPL_TS=$(echo "${TS_TPL_1}-${TS_TPL_0}"|bc -l)
    dTPL_TM=$(echo "${dTPL_TS}/60.0"|bc -l)
    dTPL_TH=$(echo "${dTPL_TS}/3600.0"|bc -l)

    dGEOS_TS=$(echo "${TS_GEOS_1}-${TS_GEOS_0}"|bc -l)
    dGEOS_TM=$(echo "${dGEOS_TS}/60.0"|bc -l)
    dGEOS_TH=$(echo "${dGEOS_TS}/3600.0"|bc -l)

    dCTEST_TS=$(echo "${TS_CTEST_1}-${TS_CTEST_0}"|bc -l)
    dCTEST_TM=$(echo "${dCTEST_TS}/60.0"|bc -l)
    dCTEST_TH=$(echo "${dCTEST_TS}/3600.0"|bc -l)

    dTS=$(echo "${TS_1}-${TS_0}"|bc -l)
    dTM=$(echo "${dTS}/60.0"|bc -l)
    dTH=$(echo "${dTS}/3600.0"|bc -l)

    printf "## CONFIG_NAME = ${CONFIG_NAME} .......................................................... \n" |& tee -a $statsf
    printf "#  Total TPL time   : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dTPL_TS $dTPL_TM $dTPL_TH |& tee -a $statsf
    printf "#  Total GEOS time  : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dGEOS_TS $dGEOS_TM $dGEOS_TH |& tee -a $statsf
    printf "#  Total CTIME time : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dCTEST_TS $dCTEST_TM $dCTEST_TH |& tee -a $statsf
    printf "#  Total build time : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dTS $dTM $dTH |& tee -a $statsf
    printf "## .................................................................................................. \n" |& tee -a $statsf
    # set +x 

fi

if [[ "${GEOS}" == "1" && "${BUILD_ONLY}" == "1" ]] || [[ "${INIT_ONLY}" == "1" ]] ; then 
    pushd $GEOSX_BUILD_DIR
    # if [ ${LCACHE} -gt 0 ]; then
    # 	pushd $croot/$PKG/$PKG/$BUILDDIR
    # else
    # 	pushd $sroot/$PKG/$PKG/$BUILDDIR
    # fi
fi

if [[ $PUBLISH_INIT == 1 ]]; then 

    echo "## $0 Note : "
    echo "#     Building $PKG modulefile ${modfile} at GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (pwd=$(pwd))"

    cat > $modfile  << EOF
#%Module1.0#####################################################################
##
## GEOS modulefile
##
## Michael E. Thomadakis, HPC R&D and Innovation, CTC
##
## Generated by $USER at installation time $(date) for
##
## PKG           : ${PKG}
## VER           : ${VER}
## HOST_CONFIG_FILENAME : ${HOST_CONFIG_FILENAME}
## ${PKG}_BRANCH   : ${GEOS_EFFECTIVE_BRANCH}
## TPL_BRANCH    : ${TPL_BRANCH}
## STACK         : ${STACK}
## ARCH          : ${ARCH}
## OS            : ${OS}
## HOST_COMPILER : ${HOST_COMPILER}
## GCC_VER       : ${GCC_VER}
## CUDA_VER      : ${CUDA_VER}
##
## Instalation location for 
## Root          : ${iroot}
## package       : ${GEOSX_DIR}}
## modulefile    : ${mroot}/${modfile}
##
## $Id: GEOS--cache-build-install.sh 2.0 2023/01/17 20:34:51 mtml Exp mtml $ ##
##
## Common modulefile part
module-whatis   "${PKG} ${VER} compiled by ${STACK} for ${ARCH} and ${OS} installed on ${GEOSX_DIR}}. "

set geosx_dir               "${GEOSX_DIR}"

setenv GEOSX_DIR            "\$geosx_dir"
setenv GEOSX_HOME           "\$geosx_dir"
setenv GEOSX_ROOT           "\$geosx_dir"
setenv HOST_CONFIG_BASE     "${HOST_CONFIG_BASE}"
setenv GEOSX_CONFIG_NAME    "${CONFIG_NAME}"
setenv HOST_CONFIG_FILENAME "${HOST_CONFIG_FILENAME}"
setenv GEOSX_VER            "${VER}"
setenv GEOS_BRANCH          "${GEOS_EFFECTIVE_BRANCH}"
setenv GEOSX_PREREQ_MODULES        "${MODULES}"
setenv GEOSX_ALL_MODULES           "${MODULES} ${mroot}/${modfile}"
setenv TPL_BRANCH                  "${TPL_EFFECTIVE_BRANCH}"
setenv GEOSX_TPL_DIR               "${GEOSX_TPL_DIR}"
setenv GEOSX_DIR                   "\${geosx_dir}"
setenv GEOSX_BUILD_DIR             "${GEOSX_BUILD_DIR}"
setenv TPL_BUILD_DIR               "${TPL_BUILD_DIR}"
setenv BUILD_TYPE		   "${BUILD_TYPE}"
setenv build_type		   "${build_type}"
setenv BUILDDIR			   "${BUILDDIR}"
setenv INSTALLDIR		   "${INSTALLDIR}"
setenv GEOSX_VENV                  "${GEOSX_VENV}"
setenv GEOSX_BUILD_VENV_DIR        "${GEOSX_BUILD_VENV_DIR}"
setenv GEOSX_INSTALL_VENV_DIR      "${GEOSX_INSTALL_VENV_DIR}"

prepend-path    DOCDIR      "\$geosx_dir"
prepend-path    DATDIR      "\$geosx_dir"

module load     ${MODULES}

prepend-path    LD_LIBRARY_PATH        "\$geosx_dir/lib"
prepend-path    PATH        "\$geosx_dir/bin"

## Auxiliary variables related to the particular GEOS/TPL build process
setenv GCC_VER              "$GCC_VER"
setenv HOST_COMPILER        "${HOST_COMPILER}"
setenv CUDA_VER             "$CUDA_VER"
setenv sroot                "$sroot"
setenv iroot                "$iroot"
setenv croot                "$croot"

setenv MPI            "$MPI"
setenv MPI_VER        "$MPI_VER"

EOF

# Specialized modulefile part for OpenMPI 
    if [ "${MPI}" == "ompi" ]; then
	cat >> $modfile <<EOF
## OpenMPI specific modulefile part

setenv OMPI_CC        $OMPI_CC  
setenv OMPI_CXX       $OMPI_CXX
setenv OMPI_F77       $OMPI_F77
setenv OMPI_FC        $OMPI_FC

# Compilers 
setenv GEOSX_CC       $GEOSX_CC
setenv GEOSX_CXX      $GEOSX_CXX
setenv GEOSX_FORT     $GEOSX_FORT

# Open MPI specifics
setenv GEOSX_MPI_DIR  $GEOSX_MPI_DIR

setenv GEOSX_MPIRUN   ${GEOSX_MPI_DIR}/bin/mpirun
setenv GEOSX_MPICC    ${GEOSX_MPI_DIR}/bin/mpicc
setenv GEOSX_MPICXX   ${GEOSX_MPI_DIR}/bin/mpic++
setenv GEOSX_MPIFORT  ${GEOSX_MPI_DIR}/bin/mpifort

# MPI independent logic
setenv MPIRUN         $GEOSX_MPIRUN

EOF
    else
	cat >> $modfile <<EOF
## IntelMPI specific modulefile part

setenv I_MPI_CC       $I_MPI_CC
setenv I_MPI_CXX      $I_MPI_CXX
setenv I_MPI_FC       $I_MPI_FC
setenv I_MPI_F90      $I_MPI_F90 


# Compilers 
setenv GEOSX_CC       $GEOSX_CC
setenv GEOSX_CXX      $GEOSX_CXX
setenv GEOSX_FORT     $GEOSX_FORT

# Open MPI specifics
setenv GEOSX_MPI_DIR  $GEOSX_MPI_DIR

setenv GEOSX_MPIRUN   ${GEOSX_MPI_DIR}/bin/mpirun
setenv GEOSX_MPICC    ${GEOSX_MPI_DIR}/bin/mpicc
setenv GEOSX_MPICXX   ${GEOSX_MPI_DIR}/bin/mpic++
setenv GEOSX_MPIFORT  ${GEOSX_MPI_DIR}/bin/mpifort

# MPI independent logic
setenv MPIRUN         $GEOSX_MPIRUN

EOF
    fi

echo "#     Publishing $PKG $modfile to ${mroot}  $mrootallu and $mrootalli"

if [ $V -ge 2 ]; then 
    set -x
fi

if [[ "${GEOS}" == "1" && "${BUILD_ONLY}" == "1" ]] || [[ "${INIT_ONLY}" == "1" ]] ; then 
    [ ! -d ${mroot} ]  && mkdir -p $mroot 
    \cp $modfile $mroot;
    # [ ! -d ${mrootallu} ]  && mkdir -p $mrootallu 
    # \cp $modfile $mrootallu
    [ ! -d ${mrootalli} ]  && mkdir -p $mrootalli 
    \cp $modfile $mrootalli
fi 

[ ! -d ./modulefiles ]  && mkdir -p ./modulefiles 
\mv $modfile ./modulefiles

if [ $V -ge 2 ]; then 
    set +x
fi



    echo "#     Building $PKG intialization file ${initrcf} "

    cat > $initrcf <<EOF
#!/bin/bash 
##
## Michael E. Thomadakis, HPC R&D and Innovation, CTC for Chevron HPC; 
##
## Initializes the base environment for $PKG $VER on ${ARCH} ${OS} built using $STACK
## by ${USER} on $(date) at $(pwd) with configuration command line options
##   [$0 "$@"]
##   Python = ${CVX_PYTHON3}
##
## PKG           : ${PKG}
## VER           : ${VER}
## ${PKG}_BRANCH   : ${GEOS_EFFECTIVE_BRANCH}
## TPL_BRANCH    : ${TPL_BRANCH}
## STACK         : ${STACK}
## ARCH          : ${ARCH}
## OS            : ${OS}
## HOST_COMPILER : ${HOST_COMPILER}
## GCC_VER       : ${GCC_VER}
## CUDA_VER      : ${CUDA_VER}
##
## GEOS built     at GEOSX_BUILD_DIR= ${GEOSX_BUILD_DIR} 
## GEOS installed at GEOSX_DIR= ${GEOSX_DIR} 
## TPL  built     at TPL_BUILD_DIR= ${GEOSX_BUILD_DIR} 
## TPL  installed at GEOSX_TPL_DIR= ${GEOSX_TPL_DIR}
## venv (built)   at GEOSX_BUILD_VENV_DIR= ${GEOSX_BUILD_VENV_DIR}
## venv (install) at GEOSX_INSTALL_VENV_DIR= ${GEOSX_INSTALL_VENV_DIR}
##
: \${V:=0}

if [ \$V -ge 4 ]; then 
    set -x
fi

if [[ "\${BASH_SOURCE[0]}" != "\${0}" ]]; then
   echo "# \${BASH_SOURCE[0]} is initializing BASH Linux environment of \$USER"
   echo "# for $PKG, $PKG $VER built with ${CONFIG_NAME} on ${ARCH} and ${OS}"  
   echo "# using s/w stack : $STACK. "  
   export INITRC=\${BASH_SOURCE[0]}
   if [ $GEOS_INSTALL -eq 0 ]; then
      echo "## Note : GEOS was built at ${GEOSX_BUILD_DIR} but not installed. "
   fi
EOF

    if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o  "${PYGEOS_BUILD}" != "0" ]; then
	echo "#     Adding support for GEOS Python and tools with virtual environment at GEOSX_INSTALL_VENV_DIR=${GEOSX_INSTALL_VENV_DIR}"

	cat >> $initrcf <<EOF
   ## Adding support for GEOS Python and tools  
   echo "# Enabling support for GEOS Python tools (python=$(which python)) with virtual environment at GEOSX_INSTALL_VENV_DIR=${GEOSX_INSTALL_VENV_DIR}"  
   unalias python python3 >& /dev/null
   export ANACONDA_VER="${ANACONDA_VER}"
   source /util/Anaconda/\${ANACONDA_VER}/etc/profile.d/conda.sh
   conda activate 
   export CVX_PYTHON3=/util/Anaconda/\${ANACONDA_VER}/bin/python
   export GEOSX_INSTALL_VENV_DIR=${GEOSX_INSTALL_VENV_DIR}
   if [ -r \$GEOSX_INSTALL_VENV_DIR/bin/activate ]; then 
   	echo "## GEOS : Activating GEOS Python virtual environment $GEOSX_VENV at installation path $GEOSX_INSTALL_VENV_DIR"
	echo "##        Using python=\$(which python) "
	source \$GEOSX_INSTALL_VENV_DIR/bin/activate
	export CVX_PYTHON3=\$(which python)
	export Python3_ROOT_DIR=\$GEOSX_INSTALL_VENV_DIR
	export Python3_EXECUTABLE=\$(basename \${CVX_PYTHON3})
	# export GEOS_CMAKE_VARS=" -DPython3_ROOT_DIR=$Python3_ROOT_DIR -DPython3_EXECUTABLE=$Python3_EXECUTABLE ${GEOS_CMAKE_VARS} "
	if [ \$V -gt 0 ]; then 
	    echo "##  Python3_ROOT_DIR = \$Python3_ROOT_DIR "
	    echo "##  Python3_EXECUTABLE = \$Python3_EXECUTABLE "
	    echo "##  CVX_PYTHON3 = \$CVX_PYTHON3 "
	    # echo "##  GEOS_CMAKE_VARS = \"\${GEOS_CMAKE_VARS}\" "
	fi
    else
	echo "## GEOS Python ERROR : Could not activate GEOS Python virtual environment $GEOSX_VENV at $GEOSX_INSTALL_VENV_DIR ! "
	echo "## Python Tools and ATS will not work ! "
	# exit 10;
    fi
    if [ "\${GEOS_INTEGRATED_TESTS}" != "0" ]; then
       echo "## GEOS Note : ATS (Integrated Tests ) can only be used from within $GEOSX_BUILD_DIR ! "
    fi		    

EOF

    fi

    cat >> $initrcf <<EOF
   ## Initializing Chevron's Environment Modules system 
   #  Selection of HPCX for CUDA 
   : \${CUDA:=$CUDA};
   export CUDA 
   source /data/saet/hpcrnd/utils/bin/modules_init.sh 
   module load "${mroot}/${modfile}"
   # module load ${MODULES} "${mroot}/${modfile}"
 else
   echo "## \${BASH_SOURCE[0]} Note : "
   echo "## to work properly, script \${BASH_SOURCE[0]} must be sourced :"
   echo "\$ source \${BASH_SOURCE[0]}"
fi

if [ \$V -ge 4 ]; then 
   set +x
fi

EOF


    chmod a+x $initrcf

    if [ $V -ge 2 ]; then 
	set -x
    fi

    echo "#     Publishing $PKG $initrcf to ${initroot} ${initrooti} "
    if [[ "${GEOS}" == "1" && "${BUILD}" == "1" ]] || [[ "${INIT_ONLY}" == "1" ]] ; then 
	[ ! -d $initroot ]  && mkdir -p $initroot
	\cp $initrcf $initroot
	# [ ! -d $initrootu ]  && mkdir -p $initrootu 
	# \cp $initrcf $initrootu
	[ ! -d $initrooti ]  && mkdir -p $initrooti
	\cp $initrcf $initrooti
    fi 
    [ ! -d .init ]  && mkdir -p ./.init
    \mv $initrcf ./.init

    if [ $V -ge 2 ]; then 
	set +x
    fi

    # INIT_ONLY != 0
    if [ ${LCACHE} -eq 3 ]; then
	echo "## GEOS : Syncing $croot/$PKG/$PKG/$BUILDDIR --> $sroot/$PKG/$PKG "
	rsync -a${RSYNCV} $croot/$PKG/$PKG/$BUILDDIR $sroot/$PKG/$PKG
    elif [  ${LCACHE} -ge 4  ]; then
	echo "## Syncing EVERYTHING ! $croot/$PKG --> $sroot "
	rsync -a${RSYNCV} $croot/$PKG $sroot
    fi

# INIT_ONLY 
fi 

#  export PYTHON3=/util/Anaconda/${CONDA_VER}/bin/python
## set(Python3_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python CACHE PATH "")
## set(Python3_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python/bin/python3 CACHE PATH "")

### gcc -march=znver4  -mdump-tune-features --help=target

### CUDA sm_80 = Ampere, sm_90 = hopper
### /usr/local/cuda/bin/nvcc
  # -gencode=arch=compute_52,code=sm_52
  # -gencode=arch=compute_60,code=sm_60
  # -gencode=arch=compute_61,code=sm_61
  # -gencode=arch=compute_70,code=sm_70
  # -gencode=arch=compute_75,code=sm_75
  # -gencode=arch=compute_80,code=sm_80
  # -gencode=arch=compute_90,code=sm_90
  # -gencode=arch=compute_90,code=compute_90
  # -O2 -o mykernel.o -c mykernel.cu


###  ## CUDA=11  sroot=${HOME}/src/GEOS_GPU WHO=CVX LCACHE=2 projdir=$(date +%F)  V=1 ccomp="gcc/11.4.0-rh8" gpu_comp="CUDA/11.8" CLONE=0 BUILD=1 TTPL_BRANCH="hypre" GGEOS_BRANCH="feature/paludettomag1/hypre-sdc" GEOS_INTEGRATED_TESTS=1 GEOS_INSTALL=0 TPL=0 GEOS=1 ./GEOS-CPU-GCC-ompi--cache-build-install.sh CPU-OPTO3-Hypre-GCC_13.2.0-ompi_hpcx-OMP

