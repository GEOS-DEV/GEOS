#
# Specify all headers
#
set(managers_headers
    DomainPartition.hpp
    EventManager.hpp
    Events/EventBase.hpp
    Events/PeriodicEvent.hpp
    Events/HaltEvent.hpp
    Events/SoloEvent.hpp
    Outputs/OutputManager.hpp
    Outputs/OutputBase.hpp
    Outputs/SiloOutput.hpp
    Outputs/RestartOutput.hpp
    Outputs/TimeHistoryOutput.hpp
    Tasks/TasksManager.hpp
    Tasks/TaskBase.hpp
    TimeHistory/TimeHistoryCollection.hpp
    TimeHistory/PackCollection.hpp
    TimeHistory/HistoryIO.hpp
    TimeHistory/HistoryDataSpec.hpp
    Outputs/BlueprintOutput.hpp
    Outputs/PythonOutput.hpp
    Functions/FunctionBase.hpp
    Functions/TableFunction.hpp
    Functions/FunctionManager.hpp
    ObjectManagerBase.hpp
    ProblemManager.hpp
    NumericalMethodsManager.hpp
    FieldSpecification/FieldSpecificationBase.hpp
    FieldSpecification/FieldSpecificationManager.hpp
    FieldSpecification/DirichletBoundaryCondition.hpp
    FieldSpecification/SourceFluxBoundaryCondition.hpp
    initialization.hpp
    GeosxState.hpp
   )


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
  list( APPEND managers_headers
        Functions/SymbolicFunction.hpp
        Functions/CompositeFunction.hpp )
  list( APPEND managers_sources
        Functions/SymbolicFunction.cpp
        Functions/CompositeFunction.cpp )
endif()

if( ENABLE_MPI )
    list( APPEND managers_headers Outputs/ChomboIO.hpp )
    list( APPEND managers_sources Outputs/ChomboIO.cpp )
endif()

if( ENABLE_VTK )
    list( APPEND managers_headers Outputs/VTKOutput.hpp )
    list( APPEND managers_sources Outputs/VTKOutput.cpp )
endif()

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

blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 HEADERS               ${managers_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX managers )
