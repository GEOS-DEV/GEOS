#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     FiniteElementDiscretizationManager.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository codingUtilities)
else()
  set( dependencyList common)
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  finiteElement
                 SOURCES               ${finiteElement_sources}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX finiteElement )

add_subdirectory( unitTests )

