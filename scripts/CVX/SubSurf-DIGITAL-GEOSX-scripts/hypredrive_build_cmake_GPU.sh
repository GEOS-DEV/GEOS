#!/bin/bash
#
# CMake GPU builds
#
# MIchael E. Thomadakis HPC INnovation and R&D 
#


if [[ ${UMPIRE} == "1" ]]; then

    cd $sroot; nowAt $0 $sroot 

    UMPIRE_DIR=${sroot}/Umpire
    UMPIRE_BUILD=${UMPIRE_DIR}/${BUILDDIR}
    UMPIRE_INSTALL=${iroot}/${ARCH}/${OS}/Umpire/${INSTALLDIR}
    # UMPIRE_INSTALL=${UMPIRE_DIR}/${INSTALLDIR}
    # UMPIRE_DIR=${CWD}/Umpire
    # UMPIRE_BUILD=${UMPIRE_DIR}/build
    # UMPIRE_INSTALL=${UMPIRE_DIR}/install-${SUFFIX}
    if [[ $INIT_ONLY != 1 ]]; then 

	echoc0 "  rm -rf ${UMPIRE_BUILD} ${UMPIRE_INSTALL}"
	rm -rf ${UMPIRE_BUILD} ${UMPIRE_INSTALL}

	if [[ $CLONE == "1" || ! -d $UMPIRE_DIR ]]; then 

	    echoc0 "  rm -rf ${UMPIRE_DIR}"
	    rm -rf ${UMPIRE_DIR}

	    echoc0 " git clone --recurse-submodules https://github.com/LLNL/Umpire.git"
	    git clone --recurse-submodules https://github.com/LLNL/Umpire.git

	fi
	# Miket: host-config cmake options
	# DCMAKE_TOOLCHAIN_FILE
	#    -DCMAKE_TOOLCHAIN_FILE=$HOST_CONFIG_FILENAME 
	echoc0 " Configuring Umpire "
	cmake \
	    -DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH \
	    -B ${UMPIRE_BUILD} -S ${UMPIRE_DIR} \
	    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
            -DUMPIRE_ENABLE_C=ON \
            -DUMPIRE_ENABLE_TOOLS=OFF \
            -DENABLE_CUDA=${ENABLE_CUDA} \
            -DENABLE_HIP=${ENABLE_HIP} \
            -DENABLE_BENCHMARKS=OFF \
            -DENABLE_EXAMPLES=OFF \
            -DENABLE_DOCS=OFF \
            -DENABLE_TESTS=OFF \
            -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
            -DCMAKE_INSTALL_LIBDIR=${UMPIRE_INSTALL}/lib \
            -DCMAKE_INSTALL_PREFIX=${UMPIRE_INSTALL}  ; RC=$? ; 
	
	cmake --build ${UMPIRE_BUILD} --parallel
	cmake --install ${UMPIRE_BUILD}
    fi
fi

# --- hypre ---
HYPRE_DIR=${sroot}/hypre
HYPRE_BUILD=${HYPRE_DIR}/${BUILDDIR}
HYPRE_INSTALL=${iroot}/${ARCH}/${OS}/hypre/${INSTALLDIR}
HYPREDRV_DIR=${sroot}/hypredrive
HYPREDRV_BUILD=${HYPREDRV_DIR}/${BUILDDIR}
HYPREDRV_INSTALL=${iroot}/${ARCH}/${OS}/hypredrive/${INSTALLDIR}

if [[ $INIT_ONLY != 1 ]]; then 
    
    if [[ $HYPRE == "1"  ]]; then 
	cd $sroot; nowAt $0 $sroot 

	echoc0 "  rm -rf ${HYPRE_BUILD} ${HYPRE_INSTALL}"
	rm -rf ${HYPRE_BUILD} ${HYPRE_INSTALL}

	if [[ $CLONE == "1" || ! -d $HYPRE_DIR ]]; then 
	    echoc0 "  rm -rf ${HYPRE_DIR} "
	    rm -rf ${HYPRE_DIR} 

	    echoc0 " git clone --depth=1 https://github.com/hypre-space/hypre.git ${HYPRE_DIR} "
	    git clone --depth=1 https://github.com/hypre-space/hypre.git ${HYPRE_DIR} || true
	fi
	# Miket: host-config cmake options
	# DCMAKE_TOOLCHAIN_FILE
	echoc0 " Configuring Hypre "
	#    -DCMAKE_TOOLCHAIN_FILE=$HOST_CONFIG_FILENAME \
	cmake \
	    -DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH \
	    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
	    -S ${HYPRE_DIR}/src -B ${HYPRE_BUILD} \
	    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	    -DCMAKE_INSTALL_PREFIX=${HYPRE_INSTALL} \
	    -DHYPRE_ENABLE_CUDA=${ENABLE_CUDA} \
	    -DHYPRE_ENABLE_HIP=${ENABLE_HIP} \
	    -Dumpire_DIR=${UMPIRE_INSTALL}/lib/cmake \
            -DMPI_C_COMPILER=${GEOSX_MPICC} \
            -DMPI_CXX_COMPILER=${GEOSX_MPICXX} \
            -DMPI_Fortran_COMPILER=${GEOSX_MPIFORT} \
            -DCMAKE_C_COMPILER=${GEOSX_CC} \
            -DCMAKE_CXX_COMPILER=${GEOSX_CXX} \
            -DCMAKE_Fortran_COMPILER=${GEOSX_FORT} \
	    -Dumpire_ROOT=${UMPIRE_INSTALL} ; RC=$? 
	    # -Dumpire_DIR=${UMPIRE_INSTALL}/lib/cmake ; RC=$? 

	    # -DCMAKE_CUDA_ARCHITECTURES=80 \
	    # -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
	    # -S ${HYPRE_DIR}/src -B ${HYPRE_BUILD} \
	    # -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	    # -DCMAKE_INSTALL_PREFIX=${HYPRE_INSTALL} \
	    # -DHYPRE_ENABLE_CUDA=${ENABLE_CUDA} \
	    # -DHYPRE_ENABLE_HIP=${ENABLE_HIP} \
	    # -Dumpire_DIR=${UMPIRE_INSTALL}/lib/cmake
	
	echoc0 " Building Hypre "
	cmake --build ${HYPRE_BUILD} --parallel
	echoc0 " Installing Hypre "
	cmake --install ${HYPRE_BUILD}
    else
	echoc0 " Did not build Hypre "
    fi
    # --- hypredrive ---
    # HYPREDRV_INSTALL=${HYPREDRV_DIR}/${INSTALLDIR}
    # HYPREDRV_DIR=${CWD}/hypredrive
    # HYPREDRV_BUILD=${HYPREDRV_DIR}/build-${SUFFIX}
    # HYPREDRV_INSTALL=${HYPREDRV_DIR}/install-${SUFFIX}

    cd $sroot; nowAt $0 $sroot 

    echoc0 "  ${HYPREDRV_DIR} ${HYPREDRV_BUILD} ${HYPREDRV_INSTALL}"
    if [[ $CLONE == "1" || ! -d $HYPREDRV_DIR ]]; then 
	echoc0 "  rm -rf ${HYPREDRV_DIR} "
	rm -rf ${HYPREDRV_DIR} 

	echoc0 " git clone --depth=1 https://github.com/hypre-space/hypredrive.git ${HYPREDRV_DIR}"
	git clone --depth=1 https://github.com/hypre-space/hypredrive.git ${HYPREDRV_DIR} || true
    fi

    echoc0 " rm -rf ${HYPREDRV_BUILD} ${HYPREDRV_INSTALL} "
    rm -rf ${HYPREDRV_BUILD} ${HYPREDRV_INSTALL}
    # cmake -S ${HYPREDRV_DIR} -B ${HYPREDRV_BUILD} 
    # Miket: host-config cmake options
    # DCMAKE_TOOLCHAIN_FILE
    echoc0 " Configuring Hypredrive "
    #    -DCMAKE_TOOLCHAIN_FILE=$HOST_CONFIG_FILENAME \
    cmake -DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH \
	-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
	-S ${HYPREDRV_DIR} -B ${HYPREDRV_BUILD} \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-DCMAKE_INSTALL_PREFIX=${HYPREDRV_INSTALL} \
	-DHYPRE_ENABLE_CUDA=${ENABLE_CUDA} \
	-DHYPRE_ENABLE_HIP=${ENABLE_HIP} \
	-Dumpire_DIR=${UMPIRE_INSTALL}/lib/cmake \
        -DMPI_C_COMPILER=${GEOSX_MPICC} \
        -DMPI_CXX_COMPILER=${GEOSX_MPICXX} \
        -DMPI_Fortran_COMPILER=${GEOSX_MPIFORT} \
        -DCMAKE_C_COMPILER=${GEOSX_CC} \
        -DCMAKE_CXX_COMPILER=${GEOSX_CXX} \
        -DCMAKE_Fortran_COMPILER=${GEOSX_FORT} \
	-DHYPRE_ROOT=${HYPRE_INSTALL} ; RC=$?;

	# -DCMAKE_CUDA_ARCHITECTURES=80 \
	# -DCMAKE_TOOLCHAIN_FILE=$HOST_CONFIG_FILENAME \
	# -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} \
	# -S ${HYPREDRV_DIR} -B ${HYPREDRV_BUILD} \
	# -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
	# -DCMAKE_INSTALL_PREFIX=${HYPREDRV_INSTALL} \
	# -DHYPRE_ROOT=${HYPRE_INSTALL} ; RC=$?;
    echoc0 " Configuring Hypredrive RC=$RC "

    echoc0 " Building Hypredrive "
    cmake --build ${HYPREDRV_BUILD} --parallel ; RC=$?;
    echoc0 " RC=$RC "
    echoc0 " Installing Hypredrive "
    cmake --install ${HYPREDRV_BUILD} ; RC=$?;
    [[ $RC == "0" ]] && PUBLISH_INIT=1
    echoc0 " RC=$RC "
   
fi

if [[ $PUBLISH_INIT == "1" || $INIT_ONLY == "1" ]]; then 

    ## Generate and install Init.rc and modulefile
    # [[ $V -ge "2" ]] && set -x ;
    set -x
    
    modfile="${ARCH}-${OS}-${PKG}-${PKG_ver}-${STACK}"; 
    initrcf="Init-${modfile}.rc";

    mrootalli="${iroot}/modulefiles/${PKG}"; 
    mroot="${HYPREDRV_INSTALL}/modulefiles";
    initrooti="${iroot}/.init/${PKG}"; 
    initroot="${HYPREDRV_INSTALL}/.init";

    # [[ $V -ge "2" ]] && set +x ;
    # set +x
    
    cd $HYPREDRV_BUILD; nowAt $0 $HYPREDRV_BUILD
    
    echoc0 "Note : "
    echoc0 "     Building $PKG modulefile ${modfile} at HYPREDRV_BUILD=$HYPREDRV_BUILD (pwd=$(pwd))"

    # set -x 
    
    cat > $modfile  <<MODFILE
#%Module1.0#####################################################################
##
## GEOS modulefile
##
## Michael E. Thomadakis, HPC R&D and Innovation, CTC
##
## Generated by $USER at installation time $(date) for
##
## PKG           : ${PKG}
## PKG_ver       : ${PKG_ver}
## HOST_CONFIG_FILENAME : ${HOST_CONFIG_FILENAME}
## STACK         : ${STACK}
## ARCH          : ${ARCH}
## OS            : ${OS}
## HOST_COMPILER : ${HOST_COMPILER}
## GCC_VER       : ${GCC_VER}
## CUDA_VER      : ${CUDA_VER}
##
## Instalation location for 
## Root          : ${iroot}
## package       : ${HYPREDRV_INSTALL}
## modulefile    : ${mroot}/${modfile}
##
##
## Common modulefile part
module-whatis   "${PKG} ${VER} compiled by ${STACK} for ${ARCH} and ${OS} installed on ${HYPREDRV_INSTALL}. "

set pkg_dir               "${HYPREDRV_INSTALL}"

setenv HYPREDRIVE_HOME             "\$pkg_dir"
setenv HYPREDRIVE_ROOT             "\$pkg_dir"
setenv HOST_CONFIG_BASE            "${HOST_CONFIG_BASE}"
setenv HOST_CONFIG_FILENAME        "${HOST_CONFIG_FILENAME}"
setenv HYPREDRIVE_VER              "${PKG_ver}"
setenv HYPREDRIVE_PREREQ_MODULES   "${GEOSX_PREREQ_MODULES}"
setenv HYPREDRIVE_ALL_MODULES      "${GEOSX_PREREQ_MODULES} ${mroot}/${modfile}"
setenv HYPREDRIVE_DIR              "${HYPREDRV_DIR}"
setenv HYPREDRIVE_BUILD_DIR        "${HYPREDRV_BUILD}"
setenv HYPRE_HOME              	   "${HYPRE_INSTALL}"
setenv HYPRE_ROOT              	   "${HYPRE_INSTALL}"
setenv HYPRE_DIR              	   "${HYPRE_DIR}"
setenv HYPRE_BUILD_DIR             "${HYPRE_BUILD}"
setenv BUILD_TYPE		   "${BUILD_TYPE}"
setenv build_type		   "${build_type}"
setenv BUILDDIR			   "${BUILDDIR}"
setenv INSTALLDIR		   "${INSTALLDIR}"

prepend-path    DOCDIR             "\$pkg_dir"
prepend-path    DATDIR             "\$pkg_dir"

module load                        ${GEOSX_PREREQ_MODULES}

# prepend-path    LD_LIBRARY_PATH        "\$pkg_dir/lib"
prepend-path    PATH               "\$pkg_dir/bin"

## Auxiliary variables related to the particular GEOS/TPL build process
setenv HOST_COMPILER        "${HOST_COMPILER}"
setenv GCC_VER              "$GCC_VER"
setenv CUDA_VER             "$CUDA_VER"
setenv sroot                "$sroot"
setenv iroot                "$iroot"
setenv croot                "$croot"

setenv MPI            "$MPI"
setenv MPI_VER        "$MPI_VER"

MODFILE

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

    if [ $V -ge 2 ]; then 
	set -x
    fi

    echo "#     Publishing $PKG $modfile to ${mroot}  $mrootallu and $mrootalli"
    [[ ! -d ${mroot} ]] && mkdir -p $mroot 
    \cp -f $modfile $mroot;
    [[ ! -d ${mrootalli} ]]  && mkdir -p $mrootalli 
    \cp -f $modfile $mrootalli

    [[ ! -d ./modulefiles ]]  && mkdir -p ./modulefiles 
    \mv -f $modfile ./modulefiles

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
## PKG_ver       : ${PKG_ver}
## STACK         : ${STACK}
## ARCH          : ${ARCH}
## OS            : ${OS}
## HOST_COMPILER : ${HOST_COMPILER}
## GCC_VER       : ${GCC_VER}
## CUDA_VER      : ${CUDA_VER}
##
## Hypredrive built at HYPREDRV_BUILD= ${HYPREDRV_BUILD} 
## Hypredrive installed at HYPREDRV_INSTALL= ${HYPREDRV_INSTALL} 
## HYPRE built at HYPRE_BUILD= ${HYPRE_BUILD} 
## HYPRE installed at HYPRE_INSTALL= ${HYPRE_INSTALL} 
##
: \${V:=0}

if [ \$V -ge 4 ]; then 
    set -x
fi

if [[ "\${BASH_SOURCE[0]}" != "\${0}" ]]; then
   echo "# \${BASH_SOURCE[0]} is initializing BASH Linux environment of \$USER"
   echo "# for $PKG version $PKG_ver built for ${ARCH} and ${OS}"  
   echo "# using s/w stack : $STACK. "  
   export INITRC=\${BASH_SOURCE[0]}

   ## Initializing Chevron's Environment Modules system 
   #  Selection of HPCX for CUDA 
   : \${CUDA:=$CUDA};
   export CUDA 
   source /data/saet/hpcrnd/utils/bin/modules_init.sh 

   module load "${mroot}/${modfile}"
 else
   echo "## ${BASH_SOURCE[0]} Note : "
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
    [[ ! -d $initroot ]]  && mkdir -p $initroot
    \cp -f $initrcf $initroot
    [[ ! -d $initrooti ]]  && mkdir -p $initrooti
    \cp -f $initrcf $initrooti
    [[ ! -d .init ]]  && mkdir -p ./.init
    \mv -f $initrcf ./.init

    if [ $V -ge 2 ]; then 
	set +x
    fi

    set +x 
    
fi


# --- Summary ---
echoc0 ""
echoc0 "=========================================================================================="
echoc0 " Build Summary"
echoc0 "=========================================================================================="
echoc0 "Build type:                ${BUILD_TYPE}"
echoc0 "hypre installation:        ${HYPRE_INSTALL}"
echoc0 "hypredrive installation:   ${HYPREDRV_INSTALL}"
if [[ -n "${UMPIRE_DIR}" ]]; then
  echoc0 "Umpire directory:          ${UMPIRE_INSTALL}"
else
  echoc0 "Umpire directory:          (not used)"
fi
echoc0 ""
echoc0 "Laplacian driver location: ${HYPREDRV_INSTALL}/bin"
# find ${HYPREDRV_INSTALL} -type f -name "laplacian*" || echoc0 "laplacian driver not found!"
echoc0 "=========================================================================================="
