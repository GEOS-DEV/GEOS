#
# Specify all headers
#
set(mainInterface_headers
    ProblemManager.hpp
    initialization.hpp
    GeosxState.hpp
   )


#
# Specify all sources
#

set(mainInterface_sources
    ProblemManager.cpp
    initialization.cpp
    GeosxState.cpp
   )

if( BUILD_OBJ_LIBS )
  set( dependencyList events dataRepository fileIO optionparser RAJA linearAlgebra conduit adiak)
else()
  set( dependencyList 
       common 
       fileIO
       linearAlgebra
       mesh
       physicsSolvers
       fileIO
       codingUtilities
       dataRepository
       finiteVolume
       finiteElement
       constitutive
       optionparser
       RAJA
       conduit )
endif()

if( ENABLE_MATHPRESSO )
   set( dependencyList ${dependencyList} mathpresso )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if ( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

if ( ENABLE_MKL )
  set( dependencyList ${dependencyList} mkl )
endif()

blt_print_target_properties( TARGET adiak )


blt_add_library( NAME                  mainInterface
                 SOURCES               ${mainInterface_sources}
                 HEADERS               ${mainInterface_headers}
                 DEPENDS_ON            adiak ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

blt_print_target_properties( TARGET mainInterface )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

#add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX mainInterface )
