set( subdirs
     math
     common
     codingUtilities
     dataRepository
     schema
     functions
     constitutive
     mesh
     linearAlgebra
     fieldSpecification
     finiteElement
     finiteVolume
     discretizationMethods
     fileIO
     physicsSolvers
     events
     mainInterface
     python
   )

if ( GEOSX_BUILD_OBJ_LIBS )
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

if( NOT GEOSX_BUILD_SHARED_LIBS AND NOT GEOSX_BUILD_OBJ_LIBS ) # if we're building full static libs and not obj_libs
  set( geosx_core_list "" )
  foreach( lib ${coreLibs} )
    list( APPEND geosx_core_list ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
  endforeach()
else( )
  set( geosx_core_list ${coreLibs} )
endif( )
set ( geosx_core_libs "${geosx_core_list}" CACHE INTERNAL "" )

set( parallelDeps "" )

if ( ENABLE_OPENMP )
  set( parallelDeps ${parallelDeps} openmp )
endif()

if ( ENABLE_CUDA )
  set( parallelDeps ${parallelDeps} cuda )
endif()

if( ENABLE_HIP )
  set( parallelDeps ${parallelDeps} blt::hip )
endif()

if( ENABLE_MPI )
  set( parallelDeps ${parallelDeps} mpi )
endif()

foreach( lib ${subdirs} )
    add_subdirectory( ${lib} )
endforeach()

if( GEOSX_BUILD_SHARED_LIBS )
  set( GEOSX_BUILD_SHARED_BOOL TRUE CACHE BOOL "" )
  set( CSL_TYPE SHARED )
else()
  set( GEOSX_BUILD_SHARED_BOOL FALSE CACHE BOOL "" )
  set( CSL_TYPE STATIC )
endif()

if( GEOSX_BUILD_OBJ_LIBS )
  message( "Using cmake  Object Libraries" )
  blt_add_library ( NAME geosx_core
                    SOURCES dummy.cpp
                    DEPENDS_ON  mainInterface physicsSolvers
                    OBJECT ${GEOSX_BUILD_OBJ_LIBS}
                    SHARED ${GEOSX_BUILD_SHARED_BOOL}
                  )
else()
  message( "Using blt blt_combine_static_libraries" )
  blt_combine_static_libraries( NAME          geosx_core
                                SOURCE_LIBS   ${coreLibs}
                                LIB_TYPE      ${CSL_TYPE}
                                LINK_PREPEND  ${GEOSX_LINK_PREPEND_FLAG}
                                LINK_POSTPEND ${GEOSX_LINK_POSTPEND_FLAG}
                              )
endif()

# Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" )
  set_target_properties( geosx_core PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS TRUE )
endif()

# To install the library with the runtime path used during the building
set_target_properties( geosx_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON )

# To install the shared library
install( TARGETS geosx_core LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX coreComponents
                       EXCLUDES cmake
                       src/coreComponents/constitutive/PVTPackage
                       src/coreComponents/mesh/PAMELA )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check 2>&1 >/dev/null"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

if(GEOSX_ENABLE_TESTS)
  add_subdirectory( unitTests )
endif( )

