#
# Specify all headers
#
set(managers_headers
    Functions/FunctionBase.hpp
    Functions/TableFunction.hpp
    Functions/FunctionManager.hpp
    FieldSpecification/FieldSpecificationBase.hpp
    FieldSpecification/FieldSpecificationManager.hpp
    FieldSpecification/DirichletBoundaryCondition.hpp
    FieldSpecification/SourceFluxBoundaryCondition.hpp
   )


#
# Specify all sources
#

set(managers_sources
    Functions/FunctionBase.cpp
    Functions/TableFunction.cpp
    Functions/FunctionManager.cpp
    FieldSpecification/FieldSpecificationBase.cpp
    FieldSpecification/FieldSpecificationManager.cpp
    FieldSpecification/DirichletBoundaryCondition.cpp
    FieldSpecification/SourceFluxBoundaryCondition.cpp
   )

if( ENABLE_MATHPRESSO )
  list( APPEND managers_headers
        Functions/SymbolicFunction.hpp
        Functions/CompositeFunction.hpp )
  list( APPEND managers_sources
        Functions/SymbolicFunction.cpp
        Functions/CompositeFunction.cpp )
endif()


if( BUILD_OBJ_LIBS )
  set( dependencyList dataRepository optionparser RAJA linearAlgebra conduit )
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

blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 HEADERS               ${managers_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX managers )
