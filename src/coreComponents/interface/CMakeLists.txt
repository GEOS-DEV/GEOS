#
# Specify all headers
#
set(interface_headers
    DomainPartition.hpp
    EventManager.hpp
    Events/EventBase.hpp
    Events/PeriodicEvent.hpp
    Events/HaltEvent.hpp
    Events/SoloEvent.hpp
    Tasks/TasksManager.hpp
    Tasks/TaskBase.hpp
    ProblemManager.hpp
    NumericalMethodsManager.hpp
    initialization.hpp
    GeosxState.hpp
   )


#
# Specify all sources
#

set(interface_sources
    DomainPartition.cpp
    EventManager.cpp
    Events/EventBase.cpp
    Events/PeriodicEvent.cpp
    Events/HaltEvent.cpp
    Events/SoloEvent.cpp
    Tasks/TaskBase.cpp
    Tasks/TasksManager.cpp
    ProblemManager.cpp
    NumericalMethodsManager.cpp
    initialization.cpp
    GeosxState.cpp
   )

if( BUILD_OBJ_LIBS )
  set( dependencyList dataRepository fileIO optionparser RAJA linearAlgebra conduit )
else()
  set( dependencyList 
       common 
       fileIO
       linearAlgebra
       mesh
       physicsSolvers
       fileIO
       mpiCommunications
       codingUtilities
       dataRepository
       finiteVolume
       meshUtilities
       finiteElement
       constitutive
       optionparser
       RAJA
       mpiCommunications
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

blt_add_library( NAME                  interface
                 SOURCES               ${interface_sources}
                 HEADERS               ${interface_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( interface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

#add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX interface )
