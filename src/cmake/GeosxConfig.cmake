set( PREPROCESSOR_DEFINES BOUNDS_CHECK
                          ADIAK
                          CALIPER
                          CHAI
                          CUDA
                          CUDA_NVTOOLSEXT
                          CUDA_STACK_SIZE
                          HIP
                          FMT_CONST_FORMATTER_WORKAROUND
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          FPE
                          HYPRE
                          LINEARALGEBRA
                          CONSTITUTIVE_DRIVERS
                          CONSTITUTIVE_MPM_ONLY
                          MATHPRESSO
                          MPM_GPU_DEBUG
                          METIS
                          MKL
                          MPI
                          PARMETIS
                          PETSC
                          PVTPackage
                          PYGEOSX
                          RAJA
                          SCOTCH
                          SEPARATION_COEFFICIENT
                          SUITESPARSE
                          SUPERLU_DIST
                          TIMERS
                          TOTALVIEW_OUTPUT
                          TRILINOS
                          VTK
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES} )
    if( ${DEP}_FOUND OR ENABLE_${DEP} OR GEOS_ENABLE_${DEP} )
        set( USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
	message(STATUS "GEOS_USE_${DEP} = ${GEOS_USE_${DEP}}")
    endif()
endforeach()

set( STRICT_PPD OPENMP )

# only activate these options if they are ENABLED AND FOUND, not if either
foreach( DEP in ${STRICT_PPD} )
    if( ${DEP}_FOUND AND ( ENABLE_${DEP} OR GEOS_ENABLE_${DEP} ) )
        set( USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
	message(STATUS "GEOS_USE_${DEP} = ${GEOS_USE_${DEP}}")
    endif()
endforeach( )

set( GEOS_USE_HYPRE_DEVICE "GEOS_USE_HYPRE_${ENABLE_HYPRE_DEVICE}" )
message( STATUS "GEOS_USE_HYPRE_DEVICE = ${GEOS_USE_HYPRE_DEVICE}")

if( GEOS_USE_CUDA AND GEOS_USE_CUDA_STACK_SIZE )
    if( CUDA_STACK_SIZE )
        if ( CUDA_STACK_SIZE MATCHES "^[0-9]+$")
            set( GEOS_CUDA_STACK_SIZE ${CUDA_STACK_SIZE} )
            message( STATUS "CUDA stack size set to ${GEOS_CUDA_STACK_SIZE}" )
        endif()
    endif()
    if( NOT GEOS_CUDA_STACK_SIZE )
        message( FATAL_ERROR "ENABLE_CUDA_STACK_SIZE activated but CUDA_STACK_SIZE not properly defined as '${CUDA_STACK_SIZE}'" )
    endif()
endif()

set( GEOS_CMAKE_BUILD_TYPE "\"${CMAKE_BUILD_TYPE}\"" )

configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp )

install( FILES ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )
