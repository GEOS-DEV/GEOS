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


echo "## $0 Note : running backend script at $(pwd) for CONFIG_NAME=$CONFIG_NAME " ;

umask 0022

: ${V:="0"} ; 
export V

if [ $V -gt 1 ]; then
    export VERBOSE=1
fi

[ "$DEBUG" == "1" -o $V -ge 10 ] && set -x

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

: ${INIT_ONLY:="0"} 
export INIT_ONLY

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
 
# LD_GEOSX_CORE : 0 = no libgeosx_core.so in LD_LIBRARY_PATH; 1 = True
: ${LD_GEOSX_CORE:="0"};
export LD_GEOSX_CORE

# GPU GEOS workaround; default : NO workaround
: ${GEOS_GPU_WORKAROUND:="0"}; 
export GEOS_GPU_WORKAROUND

# Additional variables for the TPL Cmake configure command lines -DVAR=val
: ${TPL_CMAKE_VARS:=""};
export TPL_CMAKE_VARS
 
# Additional variables for the GEOS Cmake configure command lines -DVAR=val
: ${GEOS_CMAKE_VARS:=""};
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
    : ${iroot:="/data/saet/${USER}/software"} ;
fi
: ${croot:="/dev/shm/${USER}/src"} ; 


: ${projdir:=""} 

## Bulding different branches (common or explicitly for TPL and/or GEOS)

if [ -n "$BRANCH" ]; then 
    : ${BRANCH_TPL:=${BRANCH}}
    : ${BRANCH_GEOS:=${BRANCH}}
fi 
if [ -n "${BRANCH_GEOS}" ]; then
    branch_geos=-$(echo $BRANCH_GEOS | tr ' ' '_' | tr '/' '__')
    branch="${branch_geos}"
fi
if [ -n "${BRANCH_TPL}" ]; then
    branch_tpl=--$(echo $BRANCH_TPL | tr ' ' '_' | tr '/' '__')
    branch+="${branch_tpl}"
fi


export BRANCH BRANCH_TPL BRANCH_GEOS

projdir+=$branch

if [ "${projdir}" != "" ] ; then
    projdir="-${projdir}"
fi
export projdir

# -------------------------------------------------------------------------- #
# HOST_CONFIG_DIR is the file path prefix to host-configs directory 
# Clean up non-standard CONFIG_NAME and reconstruct "proper" full file path
export HOST_CONFIG_BASE=$(basename ${CONFIG_NAME} .cmake)
export HOST_CONFIG_DIR=$(readlink -f ${HOST_CONFIG_DIR})
export HOST_CONFIG_FILENAME="${HOST_CONFIG_DIR}/${HOST_CONFIG_BASE}.cmake"

# Default environment path component names
: ${ARCH:="x86_64"} ;
: ${OS:="RHEL8"} ;
: ${PKG:="GEOS"} ;
: ${VER:="0.2.0"} ;
# : ${VER:="0.2.0${projdir}"} ;
  STACK="${HOST_CONFIG_BASE}-${build_type}"
  BUILDIR="build-${STACK}"
  INSTALLDIR="install-${STACK}"
  export GEOSX_VENV="${BUILDIR}--$GEOSX_VENV" BUILDIR INSTALLDIR
# -------------------------------------------------------------------------- #

# set umask 
umask 0022 

# -------------------------------------------------------------------------- #
# Create sroot or set it to pwd and manage source directories
# echo "## $0 PRE Note : sroot=$sroot "  ;
if [ ! -d $sroot ]; then
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

export TPL_SOURCES_DIR=${sroot}/$PKG/$TPL_PKG
export GEOSX_SOURCES_DIR=${sroot}/$PKG/$PKG

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
echo "## $0 Note : Using srcpath=$srcpath  "  ;

# Where sources are actually cloned 
export TPL_SOURCES_DIR="${sroot}/${PKG}/${TPL_PKG}"
export GEOSX_SOURCES_DIR="${sroot}/${PKG}/${PKG}"


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
if [[ "${LCACHE}" == "1" || "${LCACHE}" == "2" ]] ; then
    broot=$croot; 
else
    broot=$sroot; 
fi
export TPL_BUILD_DIR=${broot}/${PKG}/${TPL_PKG}/$BUILDIR
export GEOSX_BUILD_DIR=${broot}/${PKG}/${PKG}/$BUILDIR


# Manage installation directories
if [ "${CHAP}" == "1" ]; then
    echo "## $0 Note : /chap public installation "  ;
    export iroot="/chap"
    echo "## $0 Note : Using iroot=$iroot "  ;

    isuffix_geosx="geos/${VER}${projdir}/${ARCH}/${OS}/${INSTALLDIR}"
    prefix_geosx="${iroot}/${isuffix_geosx}" ; # [ ! -d $prefix_geosx ] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
    isuffix_tpl="TPL"
    prefix_tpl="${prefix_geosx}/${isuffix_tpl}" ; # [ ! -d $prefix_tpl ] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);
    mrootalli="${iroot}/geos/modulefiles"; # [ ! -d $mrootalli ] && mkdir -p $mrootalli; # mrootalli=$(readlink -f $mrootalli);
    initrooti="${iroot}/geos/.init"; # [ ! -d $initrooti ] && mkdir -p $initrooti; #  initrooti=$(readlink -f $initrooti);

else
    echo "## $0 Note : regular installation "  ;
    if [ ! -d $iroot ]; then
	mkdir -p $iroot ; RC=$?
	if [ $RC -ne 0 ]; then
	    echo "## $0 Warning : I cannot create $iroot ! Setting iroot = $(pwd) "  ;
	    export iroot=$(pwd) 
	fi
    fi
    # iroot=$(readlink -f $iroot) ;
    echo "## $0 Note : Using iroot=$iroot "  ;

    isuffix_tpl="${ARCH}/${OS}/${PKG}TPL/${VER}${branch_tpl}/${INSTALLDIR}"
    prefix_tpl="${iroot}/${isuffix_tpl}" ; # [ ! -d $prefix_tpl ] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);

    isuffix_geosx="${ARCH}/${OS}/${PKG}/${VER}${projdir}/${INSTALLDIR}"
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
fi

# Installation directories
# Check and create TPL installation path
if [[ "${TPL}" == "1" && "${BUILD}" == "1" ]]; then
    [ ! -d $prefix_tpl ] && mkdir -p $prefix_tpl;  # prefix_tpl=$(readlink -f $prefix_tpl);
fi

# Check and create GEOS installation path
if [[ "${GEOS}" == "1" && "${BUILD}" == "1" && "${GEOS_INSTALL}" != "0" ]]; then 
    [ ! -d $prefix_geosx ] && mkdir -p $prefix_geosx;  # prefix_geosx=$(readlink -f $prefix_geosx);
fi

# Allow to use existing TPL installation possibly from a common location by setting it to GEOSX_TPL_DIR
: ${GEOSX_TPL_DIR:=${prefix_tpl}} 
export GEOSX_TPL_DIR

# GEOSX_DIR is computed
export GEOSX_DIR=${prefix_geosx}


if [ ${V} -ge 1 ]; then 
    echo "## $0 Note :  
#    prefix_geosx=$prefix_geosx ; isuffix_geosx=$isuffix_geosx 
#    prefix_tpl=$prefix_tpl ; isuffix_tpl=$isuffix_tpl"
fi

# -------------------------------------------------------------------------- #
modfile="${ARCH}-${OS}-${PKG}-${VER}${projdir}-${STACK}"; 
# mrootallu="${HOME}/modulefiles/${PKG}"; # [ ! -d $mrootallu ] && mkdir -p $mrootallu; # mrootallu=$(readlink -f $mrootallu);

mroot="${prefix_geosx}/modulefiles"; # [ ! -d $mroot ] && mkdir -p $mroot;  # mroot=$(readlink -f $mroot);

initrcf="Init-${modfile}.rc";
initrootu="${HOME}/.init/${PKG}"; # [ ! -d $initrootu ] && mkdir -p $initrootu; # initrootu=$(readlink -f $initrootu);

# -------------------------------------------------------------------------- #
export broot croot sroot iroot mroot mrootalli initrootu initrooti isuffix_tpl isuffix_geosx prefix_geosx ARCH OS PKG VER STACK
# export broot croot sroot iroot mroot mrootallu mrootalli initrootu initrooti isuffix_tpl isuffix_geosx prefix_geosx ARCH OS PKG VER STACK


echo "# ==========================================================================================="
echo "# Automated TPL and GEOS clone-configure-build-install"
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
echo "# TPL_SOURCES_DIR      = ${TPL_SOURCES_DIR} ; where TPL sources were cloned"
echo "# GEOSX_SOURCES_DIR    = ${GEOSX_SOURCES_DIR} ; where GEOS sources were cloned"
echo "# TPL_BUILD_DIR        = ${TPL_BUILD_DIR} ; where TPL will be built"
echo "# GEOSX_BUILD_DIR      = ${GEOSX_BUILD_DIR} ; where GEOS will be built"
echo "# GEOSX_TPL_DIR        = ${GEOSX_TPL_DIR} ; if set, it will be used to let GEOS targets link against TPL at this existing location "
echo "# GEOSX_DIR            = ${GEOSX_DIR} ; where GEOS will be installed"
echo "# GEOSX_VENV           = ${GEOSX_VENV} ; name of virtual environment for Python packages"
echo "# GEOSX_INSTALL_VENV_DIR it will appear later before the virtual environment for the Python tools gets built "
echo "# CVX_PYTHON3          = ${CVX_PYTHON3} "
echo "# Scrpath dir          = ${scrpath}" 
echo "# Script dir/name      = ${scrdir}/$(basename $0) "

if [ ${V} -gt 0 ]; then
    source ${scrdir}/GEOS-quick-help.sh 
    # echo "# --------------------------------------------------------------------------------------------"
    # echo "# Script dir/name = ${scrdir}/$(basename $0) "
    # echo "# ARCH            = ${ARCH} "
    # echo "# OS              = ${OS} "
    # echo "# CLONE           = ${CLONE} : 0 = do not clone, 1 = clone "
    # echo "# BUILD_ONLY      = ${BUILD_ONLY} : 0 = clone, 1 = buld TPL or GEOS "
    # echo "# BUILD_TYPE      = ${BUILD_TYPE} ;  build_type      = ${build_type} "
    # echo "# TPL             = ${TPL} : process TPL tasks"
    # echo "# gitTPL          = ${gitTPL} : change this to point to another repo, such as Chevron's ADO"
    # echo "# TPL_UPDATE      = ${TPL_UPDATE} : 0 = do not update from Github repos; 1 = update" 
    # echo "# TPL_BUILD_ONLY  = ${TPL_BUILD_ONLY} : rebuild only these coma-separated TPL packages without re-configuring TPL (0: build all TPL)"
    # echo "# TPL_RESET       = ${TPL_RESET} : '--type { HEAD~n | hash }' to reset before building TPL targets'" 
    # echo "# TPL_CMAKE_VARS  = \"${GEOS_CMAKE_VARS}\" ; additional TPL CMake configure variables" 
    # echo "# GEOS            = ${GEOS} : process GEOS tasks" 
    # echo "# gitGEOS         = ${gitGEOS} : change this to point to another repo, such as Chevron's ADO"
    # echo "# GEOS_UPDATE     = ${GEOS_UPDATE} : 0 = do not update from Github repos; 1 = update" 
    # echo "# GEOS_REBUILD    = ${GEOS_REBUILD} : 0 = configure and build all GEOS targets, 1 = just rebuild" 
    # echo "# GEOS_BUILD_ONLY = ${GEOS_BUILD_ONLY} : 0 = build all targets ; 1 = only build geosx ; 'list of build targets'" 
    # echo "# GEOS_RESET      = ${GEOS_RESET} : '--type { HEAD~n | hash }' to reset before building GEOS targets'" 
    # echo "# GEOS_CMAKE_VARS = \"${GEOS_CMAKE_VARS}\" ; additional GEOS CMake configure variables" 
    # echo "# GEOS_MAIN_BRANCH = ${GEOS_MAIN_BRANCH} : check out this GEOS branch before processing GEOS tasks" 
    # echo "# BRANCH_TPL       = ${BRANCH_TPL} : check out this TPL branch before processing TPL tasks" 
    # echo "# BRANCH_GEOS      = ${BRANCH_GEOS} : check out this GEOS branch before processing GEOS tasks" 
    # echo "# --------------------------------------------------------------------------------------------"
    # env | egrep 'ENABLE(.*)=' | sort | gawk '{print "# ", $0 }'
    # echo "# ENABLE_CALIPER_HYPRE    = ${ENABLE_CALIPER_HYPRE} : enable detailed Hypre profiling "
    # echo "# ENABLE_CUDA_NVTOOLSEXT  = ${ENABLE_CUDA_NVTOOLSEXT} : enable CUDA performance tracing "
    echo "# --------------------------------------------------------------------------------------------"
    env | egrep -i 'slurm(.*)=' | sort | gawk '{print "# ", $0 }'
    echo "# --------------------------------------------------------------------------------------------"
    echo "# prefix_tpl      = ${prefix_tpl} "
    echo "# prefix_geosx    = ${prefix_geosx} "
    echo "# projdir         = ${projdir}" 
    echo "# modfile         = ${modfile}" 
    echo "# sroot           = ${sroot}" 
    echo "# croot           = ${croot}" 
    echo "# broot           = ${broot}" 
    echo "# iroot           = ${iroot}" 
    echo "# srcpath         = ${srcpath}" 
    echo "# BUILDIR         = ${BUILDIR}" 
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
pushd $srcpath
echo "## Now at $(pwd) ( must == $srcpath )"

if [ "${INIT_ONLY}" != "1" ]; then 

    TS_0=$(date "+ %s.%N")

    if [ "${GEOS_REBUILD}" != "1"   ] ; then 

	if [ "${CLONE}" == "1" ]; then
	    echo "## GEOS : Cloning $gitGEOS --> $(pwd)/$PKG "
	    git clone ${gitGEOS} $PKG
	fi

	if [ "${GEOS_UPDATE}" == "1"  ] ; then 
	    pushd ${PKG}
	    echo "## GEOS : Updating ${PKG} at $(pwd) (also update the GEOS/host-configs before configuring and building TPL) "
	    echo "## Origin : ${gitGEOS}"

	    git checkout $GEOS_MAIN_BRANCH

	    git-lfs install

	    git pull $GEOS_PULL_OPTIONS

	    git lfs pull

	    git submodule init

	    # Building integratedTests now mandatory
	    # if [ "$GEOS_INTEGRATED_TESTS" == "0" ] ; then 
		# git submodule deinit integratedTests
	    # fi

	    git submodule update

	    ## Detailed submodule updates 
	    git submodule update --init src/coreComponents/LvArray
	    git submodule update --init src/cmake/blt
	    git submodule update --init src/coreComponents/constitutive/PVTPackage
	    git submodule update --init src/coreComponents/fileIO/coupling/hdf5_interface
	    # git submodule update --init src/externalComponents/PAMELA

	    # Populate $srcpath/$PKG/host-configs/CVX with CVX specific host_configs 
	    echo "## GEOS : Populating CVX specific $scrdir/$PKG/host-configs/CVX  --> $srcpath/$PKG/host-configs/  "
	    rsync -a${RSYNCV} $scrdir/$PKG/host-configs/CVX $srcpath/$PKG/host-configs/
	    # echo "## Copying $srcpath/$PKG/host-configs/tpls.cmake --> $srcpath/$PKG/host-configs/CVX/tpls.cmake  "
	    # \cp -f $srcpath/$PKG/host-configs/tpls.cmake $srcpath/$PKG/host-configs/CVX/tpls.cmake 
	    ls -l $srcpath/$PKG/host-configs/tpls.cmake $srcpath/$PKG/host-configs/CVX/tpls.cmake | gawk '{print "#   ", $0}'
	    echo "## GEOS : Updating local host-configs with latest public GihHub : $srcpath/$PKG/host-configs --> $scrdir/$PKG/ excluding our's at ./CVX/* "
	    rsync -a${RSYNCV} $srcpath/$PKG/host-configs --exclude '*CVX/*' $scrdir/$PKG/ 
	fi
    fi

    ## TPL clone, configure and build-install tasks -------------------------------------------------------------------------
    # pushd $srcpath
    echo "## TPL : Now at $(pwd) ( must == $srcpath )"

    if [ "${TPL}" == "1" ] ; then 

	TS_TPL_0=$(date "+ %s.%N")
	
	echo "## TPL : Bringing in ${gitTPL} --> ${TPL_PKG}"

	if [ "${CLONE}" == "1" ]; then
	    git clone ${gitTPL} ${TPL_PKG}
	    # https://github.com/GEOS/thirdPartyLibs.git
	fi

	cd ${TPL_PKG}
	if [ "${TPL_UPDATE}" == "1"   ] ; then 
	    # Update the TPL packages 
	    echo "## TPL : Updating ${TPL_PKG} at $(pwd) "
	    echo "## Origin : ${gitTPL}"

	    git-lfs install

	    git pull $TPL_PULL_OPTIONS

	    git lfs pull

	    git submodule init
	    git submodule update
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

		echo "## TPL : Syncing $scrdir/$PKG/$PKG/host-configs --> $croot/$PKG/$PKG "
		rsync -a${RSYNCV} $scrdir/$PKG/$PKG/host-configs $croot/$PKG/$PKG

		pushd $croot/$PKG/${TPL_PKG} 
		echo "## Now at $(pwd) "
	    fi

	    # BUILDIR="build-${CONFIG_NAME}-${build_type}"
	    if [ ! -d ${BUILDIR} ]; then
		echo "## TPL : Creating ${BUILDIR} at $(pwd) "
		mkdir -p ${BUILDIR}; 
	    fi

	    ## Checkout required TPL branch to work with
	    if [ -n "${BRANCH_TPL}" ]; then
		
		echo "## TPL : Checking out BRANCH=${BRANCH_TPL} "
		git checkout $BRANCH_TPL ; RC=$?
		if [ $RC -ne 0 ]; then
		    echo "## $0 TPL Error : There is no TPL branch $BRANCH_TPL at $(pwd) ! Exititng ... "
		    exit $RC;
		fi
		echo "## TPL : git pull BRANCH=${BRANCH_TPL} "
		git pull
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
		    -bt ${BUILD_TYPE} -DGEOSX_TPL_DIR=$GEOSX_TPL_DIR -ip $GEOSX_TPL_DIR -DNUM_PROC=${NUM_PROC} |& tee Config-Build-${DT_GEOSX}.log  
 

		pushd ./${BUILDIR}
		export TPL_BUILDIR_PATH=$(readlink -f $(pwd))
		echo "## TPL : Building thirdPartyLibs at TPL_BUILDIR_PATH=$TPL_BUILDIR_PATH (pwd=$(pwd)) "
		make -k |& tee Make-${DT_GEOSX}.log 
	    else 
		# Only rebuild list of TPL targets
		pushd ./${BUILDIR}
		echo "## TPL : Re-building TPL packages without reconfiguration at $(pwd) (TPL_BUILD_ONLY=${TPL_BUILD_ONLY})"
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
$(env | egrep 'fix|Python3|PYTH|ENV|ENABLE(.*)|GEOS|CONF|CACHE|CUDA|ARCH|OS|PKG|VER|STACK|HPCX|MPI|STACK' | sort | gawk '{print "# ", $0 }')" |& tee -a $GEOSX_TPL_DIR/${build_configuration}
	fi

	TS_TPL_1=$(date "+ %s.%N")

    fi 

    ## GEOS clone, configure, build and install tasks  ---------------------------------------------------------------------------------
    if [ "${GEOS}" == "1" ] ; then 
	TS_GEOS_0=$(date "+ %s.%N")

	if [ "${BUILD_ONLY}" == "1" ]; then
	    ## Configure and Build 
	    echo "## GEOS : Building ${PKG}"

	    pushd $sroot/${PKG}/${PKG}
	    echo "## Now at $(pwd) "

	    ## Checkout required GEOS branch to work with
	    if [ -n "${BRANCH_GEOS}" ]; then

		echo "## GEOS : Checking out BRANCH=${BRANCH_GEOS} "
		git checkout $BRANCH_GEOS ; RC=$?
		if [ $RC -ne 0 ]; then
		    echo "## $0 GEOS Error : There is no GEOS branch $BRANCH_GEOS at $(pwd) ! Exititng ... "
		    exit $RC;
		fi
		if [ "${GEOS_UPDATE}" == "1" ] ; then 
		    echo "## GEOS : git pull BRANCH=${BRANCH_GEOS} "
		    git pull
		    git-lfs pull
		fi
	    fi

	    if [ -n "${GEOS_RESET}" ]; then 
		echo "## GEOS : git reset ${GEOS_RESET} "
		git reset ${GEOS_RESET} 
	    fi

	    if [ ${LCACHE} -eq 2 ]; then
		echo "## GEOS : Clearing $croot/$PKG/$PKG "
		\rm -fr $croot/$PKG/$PKG
	    fi

	    if [ ${LCACHE} -gt 0 ]; then
		echo "## GEOS : Syncing $sroot/$PKG/$PKG --> $croot/$PKG "
		if [[ "$GEOS_INTEGRATED_TESTS" == "1"  ]]; then 
		    rsync --exclude "$PKG/*build*" --exclude "$PKG/*install*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		else
		    rsync --exclude "$PKG/*build*" --exclude "$PKG/*install*" --exclude "$PKG/*inputFiles*"  --exclude "$PKG/*integratedTests*" -a${RSYNCV} $sroot/$PKG/$PKG  $croot/$PKG
		fi
		pushd $croot/$PKG/$PKG 
	    fi

	    ## GEOS Configuration step
	    if [ ! -d ${BUILDIR} ]; then
		echo "## GEOS : Creating BUILDIR=${BUILDIR} at $(pwd) "
		mkdir -p ${BUILDIR}; 
	    fi
	    echo "## GEOS : BUILDIR=${BUILDIR} at $(pwd) "
	    ## MikeT: Check if the GEOSX_BUILD_DIR is the same as the one already set
	    echo "## GEOS : GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (Pre)"
	    export GEOSX_BUILD_DIR=$(readlink -f ${BUILDIR})
 	    echo "## GEOS : GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (Post)"

	    # GEOS_REBUILD == 1 does not install/rebuild ATS or Python Tools
	    if [ "${GEOS_REBUILD}" != "1"   ] ; then 

		if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o "${PYGEOS_BUILD}" != "0" ]; then
		    # Prepare virtual environment for Python tools and Designate compiler for mpi4py 
		    # GEOSX_BUILD_VENV

		    # Setup virtual environment as ../$BUILDIR--$GEOSX_VENV a peer dir to $BUILDIR
		    export GEOSX_BUILD_VENV_DIR=${GEOSX_TPL_DIR}/$GEOSX_VENV
		    export GEOSX_INSTALL_VENV_DIR=${GEOSX_BUILD_VENV_DIR}
		    
		    echo "## GEOS : Preparing GEOS Python virtual environment $GEOSX_VENV  (at $(pwd)) [GEOSX_BUILD_VENV_DIR=$GEOSX_BUILD_VENV_DIR]"
		    echo "##        Installing GEOS Python virtual environment at $GEOSX_BUILD_VENV_DIR [GEOSX_INSTALL_VENV_DIR=$GEOSX_INSTALL_VENV_DIR]"
		    # $CVX_PYTHON3 -m pip install venv
		    $CVX_PYTHON3 -m venv $GEOSX_BUILD_VENV_DIR
		    if [ -r $GEOSX_BUILD_VENV_DIR/bin/activate ]; then 
			source $GEOSX_BUILD_VENV_DIR/bin/activate
			RC=$?
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
		
		echo "## GEOS : Configuring GEOS at $(pwd) using Python=$CVX_PYTHON3"
		
		$ECHO $CVX_PYTHON3 scripts/config-build.py "$GEOS_CMAKE_VARS" -hc $HOST_CONFIG_FILENAME -bt ${BUILD_TYPE} -ip $GEOSX_DIR -DGEOSX_DIR=$GEOSX_DIR -DGEOSX_TPL_DIR=$GEOSX_TPL_DIR |& tee Config-${DT_GEOSX}.log 

	    else
		echo "## GEOS Note : At this state a 'configure' step must have already taken place for ${CONFIG_NAME} "
	    fi
	    
	    if [ "${GEOS_BUILD_ONLY}" == "1" ]; then
		export GEOS_BUILD_ONLY="geosx"
	    elif [ "${GEOS_BUILD_ONLY}" == "0" ]; then
		export GEOS_BUILD_ONLY="all"
	    fi

	    pushd ./${BUILDIR}

	    ## GEOS Build step
	    echo "## GEOS : Building GEOS targets ${GEOS_BUILD_ONLY} at $(pwd) "
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
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make geosx_python_tools_clean
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make geosx_python_tools
	    fi

	    if [ "${GEOS_INTEGRATED_TESTS}" != "0" ]; then
		# export PATH="${HOME}/.local/bin:${PATH}"
		echo "## GEOS : Building intregratedTests (ats_environment) at $GEOSX_BUILD_VENV_DIR "
		# export LD_LIBRARY_PATH="/usr/lib64:/usr/lib:${LD_LIBRARY_PATH}"
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make ats_environment
		VERBOSE=1 env  MPICC=$GEOSX_MPICC make ats_clean
		###  OMP_NUM_THREADS=1 ./geos_ats.sh --machine openmpi --ats openmpi_args=--oversubscribe --ats openmpi_procspernode=4 --ats openmpi_maxprocs=8 --ats openmpi_args=--report-bindings --ats openmpi_args="--bind-to none" --ats openmpi_mpirun $(which mpirun) --ats openmpi_install $HPCX_MPI_DIR
	    fi

	    ## Installation Phase
	    if [ "${GEOS_INSTALL}" == "1" ]; then 
		echo "## GEOS : Installing GEOS in $GEOSX_DIR "
		ln -s $GEOSX_DIR
		make -k install |& tee -a Make-install-${DT_GEOSX}.log 
		cp $GEOSX_TPL_DIR/${build_configuration} $GEOSX_DIR

		## Check for installation issues
		if [ ! -x $GEOSX_DIR/bin/geosx -o "${GEOS_GPU_WORKAROUND}" == "1" ]; then
		    echo "## GEOS Note : Install workaround "    |& tee Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./include $GEOSX_DIR  |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./bin $GEOSX_DIR      |& tee -a Make-install-${DT_GEOSX}.log 
		    cp  ./bin/geosx $GEOSX_DIR/bin/geos     |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./tests/* $GEOSX_DIR/bin |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./lib $GEOSX_DIR      |& tee -a Make-install-${DT_GEOSX}.log 
		    rsync -a${RSYNCV} ./share $GEOSX_DIR    |& tee -a Make-install-${DT_GEOSX}.log 
		    mkdir -p $GEOSX_DIR/lib64               |& tee -a Make-install-${DT_GEOSX}.log 
		    cp ./lib/libbenchmark* $GEOSX_DIR/lib64 |& tee -a Make-install-${DT_GEOSX}.log 
		    $GEOSX_DIR/bin/geosx -s $GEOSX_DIR/this-schema.xsd |& tee -a Make-install-${DT_GEOSX}.log 
		fi
		# Propagate Python's virtual environment properly to install location
		# ATS may only be used at the build location !
		# if [ ${GEOS_PYTHON_TOOLS} -gt 0 -o "${GEOS_INTEGRATED_TESTS}" != "0" -o  "${PYGEOS_BUILD}" != "0" ]; then
		#     echo "## GEOS : Installing Python Tools GEOS in $GEOSX_INSTALL_VENV_DIR "
		#     echo "## GEOS : Copying python virtual environment $GEOSX_VENV at $GEOSX_BUILD_VENV_DIR to $GEOSX_INSTALL_VENV_DIR "
		#     rsync -a${RSYNCV} ${GEOSX_BUILD_VENV_DIR} $GEOSX_DIR/..
		#     echo "## GEOS NOTE : ATS may only be used at the build location in $GEOSX_BUILD_VENV_DIR !"
		#     echo "##             Other Python tools could be used out of the install location in GEOSX_INSTALL_VENV_DIR !"
		# fi
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
$(env | egrep 'fix|Python3|PYTH|ENV|ENABLE(.*)|GEOS|CONF|CACHE|CUDA|ARCH|OS|PKG|VER|STACK|HPCX|MPI|STACK' | sort | gawk '{print "# ", $0 }') " |& tee -a $GEOSX_DIR/${build_configuration}
	fi

	TS_GEOS_1=$(date "+ %s.%N")

	if [ ${LCACHE} -ge 3 ]; then
	    echo "## GEOS : Syncing $croot/$PKG/$PKG/$BUILDIR --> $sroot/$PKG/$PKG "
	    rsync -a${RSYNCV} $croot/$PKG/$PKG/$BUILDIR $sroot/$PKG/$PKG
	fi

    fi 

    ## Run ctest -V  

    if [ ${CTEST} -gt 0 ]; then
	TS_CTEST_0=$(date "+ %s.%N")
	pushd $GEOSX_BUILD_DIR

	# if [ ${LCACHE} -gt 0 ]; then
	#     pushd $croot/$PKG/$PKG/$BUILDIR
	# else
	#     pushd $sroot/$PKG/$PKG/$BUILDIR
	# fi
	echo "## GEOS : Running ctest -V on $BUILDIR at $GEOSX_BUILD_DIR. If there is no $BUILDIR, this will fail."
	ctest -V |& tee CTEST-${CONFIG_NAME}_${DT_GEOSX}.log 

	echo "## GEOS : Syncing $croot/$PKG/$PKG/$BUILDIR/Testing_${DT_GEOSX} --> $sroot/$PKG/$PKG/$BUILDIR"
	if [ ${LCACHE} -gt 0 ]; then
	    mkdir -p $sroot/$PKG/$PKG/$BUILDIR
	    mv $croot/$PKG/$PKG/$BUILDIR/Testing $croot/$PKG/$PKG/$BUILDIR/Testing_${DT_GEOSX}
	    rsync -a${RSYNCV} $croot/$PKG/$PKG/$BUILDIR/Testing_${DT_GEOSX} $sroot/$PKG/$PKG/$BUILDIR    
	else
	    mv $sroot/$PKG/$PKG/$BUILDIR/Testing $sroot/$PKG/$PKG/$BUILDIR/Testing_${DT_GEOSX}
	fi 
	TS_CTEST_1=$(date "+ %s.%N")
    fi
    
    echo "## ${CONFIG_NAME} End  : $(date +%F-%H%M%S) $(date +%s.%N) "

    # 

    if [  ${LCACHE} -ge 4  ]; then
	echo "## Syncing EVERYTHING ! $croot/$PKG --> $sroot "
	rsync -a${RSYNCV} $croot/$PKG $sroot
    fi

    # set -x 
    TS_1=$(date "+ %s.%N")

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

    printf "## CONFIG_NAME = ${CONFIG_NAME} .......................................................... \n" 
    printf "## Total TPL time   : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dTPL_TS $dTPL_TM $dTPL_TH
    printf "## Total GEOS time  : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dGEOS_TS $dGEOS_TM $dGEOS_TH
    printf "## Total CTIME time : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dCTEST_TS $dCTEST_TM $dCTEST_TH 
    printf "## Total build time : %12.3lf sec;  %12.3lf mins; %12.3lf hours \n" $dTS $dTM $dTH 
    printf "## .................................................................................................. \n" 
    # set +x 

fi

if [[ "${GEOS}" == "1" && "${BUILD_ONLY}" == "1" ]] || [[ "${INIT_ONLY}" == "1" ]] ; then 
    pushd $GEOSX_BUILD_DIR
    # if [ ${LCACHE} -gt 0 ]; then
    # 	pushd $croot/$PKG/$PKG/$BUILDIR
    # else
    # 	pushd $sroot/$PKG/$PKG/$BUILDIR
    # fi
fi

echo "## $0 Note : building $PKG modulefile ${modfile} at $GEOSX_BUILD_DIR=$GEOSX_BUILD_DIR (pwd=$(pwd))"

cat > $modfile  << EOF
#%Module1.0#####################################################################
##
## GEOS modulefile
##
## Michael E. Thomadakis, HPC R&D and Innovation, CTC
##
## Generated at installation time on $(date) for $PKG $VER on ${ARCH} ${OS} using $STACK by $USER
##
## $Id: GEOS--cache-build-install.sh 2.0 2023/01/17 20:34:51 mtml Exp mtml $ ##
##
## Instalation location for 
## Root:       ${iroot}
## package:    ${GEOSX_DIR}}
## modulefile: ${mroot}/${modfile}
##

## Common modulefile part
module-whatis   "${PKG} ${VER} compiled by ${STACK} for ${ARCH} and ${OS} installed on ${GEOSX_DIR}}. "

set geosx_dir               "${GEOSX_DIR}"

setenv GEOSX_HOME           "\$geosx_dir"
setenv GEOSX_ROOT           "\$geosx_dir"
setenv GEOSX_CONFIG_NAME    ${CONFIG_NAME}
setenv GEOSX_VER            ${VER}
setenv GEOSX_PREREQ_MODULES        "${MODULES}"
setenv GEOSX_ALL_MODULES           "${MODULES} ${mroot}/${modfile}"
setenv GEOSX_TPL_DIR               "${GEOSX_TPL_DIR}"
setenv GEOSX_DIR                   "\${geosx_dir}"
setenv GEOSX_BUILD_DIR             "${GEOSX_BUILD_DIR}"
setenv TPL_BUILD_DIR               "${TPL_BUILD_DIR}"
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
setenv sroot                $sroot
setenv iroot                $iroot
setenv croot                $croot

setenv MPI            $MPI
setenv MPI_VER        $MPI_VER

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

echo "## $0 Note : Publishing $PKG $modfile to ${mroot}  $mrootallu and $mrootalli"

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



echo "## $0 Note : building $PKG intialization file ${initrcf} "

cat > $initrcf <<EOF
#!/bin/bash 
##
## Michael E. Thomadakis, HPC R&D and Innovation, CTC for Chevron HPC; 
##
## Initializes the base environment for $PKG $VER on ${ARCH} ${OS} built using $STACK
## by ${USER} on $(date) at $(pwd) with configuration command line options
##   [$0 "$@"]
##   Python = ${CVX_PYTHON3}
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
    echo "## Adding support for GEOS Python and tools with virtual environment at GEOSX_INSTALL_VENV_DIR=${GEOSX_INSTALL_VENV_DIR}"

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
   : \${CUDA:=$CUDA};
   export CUDA 
   source /data/saet/hpcrnd/utils/bin/modules_init.sh 
   module load "${mroot}/${modfile}"
   # module load ${MODULES} "${mroot}/${modfile}"
 else
   echo "## $0 Note : "
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
# set -x

echo "## Note : Publishing $PKG $initrcf to ${initroot} ${initrootu} ${initrooti} "
if [[ "${GEOS}" == "1" && "${BUILD_ONLY}" == "1" ]] || [[ "${INIT_ONLY}" == "1" ]] ; then 
    [ ! -d $initroot ]  && mkdir -p $initroot
    \cp $initrcf $initroot
    [ ! -d $initrootu ]  && mkdir -p $initrootu  
    \cp $initrcf $initrootu
    [ ! -d $initrooti ]  && mkdir -p $initrooti
    \cp $initrcf $initrooti
fi 
[ ! -d .init ]  && mkdir -p ./.init
\mv $initrcf ./.init

if [ $V -ge 2 ]; then 
    set +x
fi

  # INIT_ONLY != 0


#  export PYTHON3=/util/Anaconda/${CONDA_VER}/bin/python
## set(Python3_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python CACHE PATH "")
## set(Python3_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python/bin/python3 CACHE PATH "")


