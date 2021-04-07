#
# Specify all sources
#

set(managers_sources
    DomainPartition.cpp
    EventManager.cpp
    Events/EventBase.cpp
    Events/PeriodicEvent.cpp
    Events/HaltEvent.cpp
    Events/SoloEvent.cpp
    Outputs/OutputManager.cpp
    Outputs/OutputBase.cpp
    Outputs/SiloOutput.cpp
    Outputs/RestartOutput.cpp
    Outputs/TimeHistoryOutput.cpp
    Outputs/BlueprintOutput.cpp
    Tasks/TaskBase.cpp
    Tasks/TasksManager.cpp
    TimeHistory/PackCollection.cpp
    Outputs/PythonOutput.cpp
    Functions/FunctionBase.cpp
    Functions/TableFunction.cpp
    Functions/FunctionManager.cpp
    ObjectManagerBase.cpp
    ProblemManager.cpp
    NumericalMethodsManager.cpp
    FieldSpecification/FieldSpecificationBase.cpp
    FieldSpecification/FieldSpecificationManager.cpp
    FieldSpecification/DirichletBoundaryCondition.cpp
    FieldSpecification/SourceFluxBoundaryCondition.cpp
    initialization.cpp
    GeosxState.cpp
   )

if( ENABLE_MATHPRESSO )
  list( APPEND managers_sources
        Functions/SymbolicFunction.cpp
        Functions/CompositeFunction.cpp )
endif()

if( ENABLE_MPI )
    list( APPEND managers_sources Outputs/ChomboIO.cpp )
endif()

if( ENABLE_VTK )
    list( APPEND managers_sources Outputs/VTKOutput.cpp )
endif()

if( BUILD_OBJ_LIBS )
  set( dependencyList dataRepository fileIO optionparser RAJA linearAlgebra conduit )
else()
  set( dependencyList common fileIO linearAlgebra optionparser RAJA mpiCommunications conduit )
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

blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX managers )
