#
# Specify all headers
#
set(events_headers
    EventManager.hpp
    EventBase.hpp
    PeriodicEvent.hpp
    HaltEvent.hpp
    SoloEvent.hpp
    tasks/TasksManager.hpp
    tasks/TaskBase.hpp
   )


#
# Specify all sources
#

set(events_sources
    EventManager.cpp
    EventBase.cpp
    PeriodicEvent.cpp
    HaltEvent.cpp
    SoloEvent.cpp
    tasks/TaskBase.cpp
    tasks/TasksManager.cpp
   )

set( dependencyList dataRepository fileIO lvarray common conduit)

if( ENABLE_MPI )
  set( dependencyList ${dependencyList} mpi )
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

blt_add_library( NAME                  events
                 SOURCES               ${events_sources}
                 HEADERS               ${events_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( events PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

#add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX events )
