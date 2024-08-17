set( subdirs
     math
     common
     codingUtilities
     dataRepository
     schema
     functions
     constitutive
     constitutiveDrivers
     mesh
     denseLinearAlgebra
     linearAlgebra
     fieldSpecification
     finiteElement
     finiteVolume
     discretizationMethods
     fileIO
     physicsSolvers
     events
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

foreach( lib ${subdirs} )
  add_subdirectory( ${lib} )
endforeach()

add_library(geosx_core INTERFACE)
if( GEOS_BUILD_SHARED_LIBS)
  target_link_libraries(geosx_core INTERFACE mainInterface )
else()
  get_target_property( libs geosx_core INTERFACE_LINK_LIBRARIES )

  set( linkLibs "" )
  foreach( lib ${coreLibs} )
    list( APPEND linkLibs "$<LINK_LIBRARY:WHOLE_ARCHIVE,${lib}>" )
  endforeach()
  message( "linkLibs = ${linkLibs}" )
  target_link_libraries(geosx_core INTERFACE ${linkLibs} )
endif()
blt_print_target_properties( TARGET geosx_core )

geosx_add_code_checks( PREFIX coreComponents
                       EXCLUDES cmake constitutive/PVTPackage )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check 2>&1 >/dev/null"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

