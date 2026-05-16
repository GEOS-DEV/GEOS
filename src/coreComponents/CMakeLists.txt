set( subdirs
     math
     common
     codingUtilities
     dataRepository )

if( GEOS_ENABLE_DENSE_LINEARALGEBRA )
  list( APPEND subdirs denseLinearAlgebra )
endif()

if( GEOS_ENABLE_LINEARALGEBRA )
  list( APPEND subdirs linearAlgebra )
endif()

list( APPEND subdirs
     functions
     constitutive )

if( GEOS_ENABLE_CONSTITUTIVE_DRIVERS )
  list( APPEND subdirs constitutiveDrivers )
endif()

list( APPEND subdirs
     schema
     finiteElement
     mesh
     fieldSpecification
     finiteVolume
     discretizationMethods
     events
     fileIO
     physicsSolvers
     mainInterface )

unset( parallelDeps )

if ( ENABLE_OPENMP )
  list( APPEND parallelDeps openmp )
endif()

if ( ENABLE_CUDA )
  list( APPEND parallelDeps cuda )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND parallelDeps CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  list( APPEND parallelDeps blt::hip )
endif()

if( ENABLE_MPI )
  list( APPEND parallelDeps mpi )
endif()

if ( GEOS_BUILD_OBJ_LIBS )
  set( LVARRAY_BUILD_OBJ_LIBS TRUE CACHE BOOL "" FORCE )
endif()

add_subdirectory( LvArray )

if( ENABLE_PYGEOSX )
  list( APPEND coreLibs pylvarray )
endif()

foreach( lib ${subdirs} )
  if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
    list( APPEND coreLibs ${lib} )
  endif()
endforeach()

foreach( subdir ${subdirs} )
  add_subdirectory( ${subdir} )
endforeach()

add_library(geosx_core INTERFACE)

if( GEOS_BUILD_SHARED_LIBS)
  target_link_libraries(geosx_core INTERFACE mainInterface )
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(geosx_core INTERFACE -Wl,--no-as-needed)
  endif()
else()
  target_link_libraries(geosx_core INTERFACE $<LINK_LIBRARY:WHOLE_ARCHIVE,mainInterface> )
endif()

geosx_add_code_checks( PREFIX coreComponents
                       EXCLUDES cmake constitutive/PVTPackage constitutive/HPCReact )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check "
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

if( GEOS_ENABLE_TESTS )
  add_subdirectory( integrationTests )
endif( )
