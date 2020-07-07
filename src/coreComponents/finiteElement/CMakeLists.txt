#
# Specify all headers
#
set( finiteElement_headers
     FiniteElementDiscretization.hpp
     ElementLibrary/FiniteElement.h
     ElementLibrary/FiniteElementBase.h
     ElementLibrary/TetrahedralElement.hpp
     Kinematics.h
     basis/BasisBase.hpp
     basis/LagrangeBasis.hpp
     basis/Polynomial.hpp
     quadrature/GaussQuadrature.hpp
     quadrature/QuadratureBase.hpp
     quadrature/QuadratureRuleManager.hpp
     kernelInterface/KernelBase.hpp
     kernelInterface/ImplicitKernelBase.hpp
     FiniteElementDiscretizationManager.hpp
     basis/BasisFunctionManager.hpp
     FiniteElementShapeFunctionKernelBase.hpp
     LinearTetrahedronShapeFunctionKernel.hpp
     BiLinearWedgeShapeFunctionKernel.hpp
     TrilinearHexahedronShapeFunctionKernel.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     ElementLibrary/FiniteElement.cpp
     ElementLibrary/FiniteElementBase.cpp
     ElementLibrary/TetrahedralElement.cpp
     basis/BasisBase.cpp
     basis/LagrangeBasis.cpp
     basis/Polynomial.cpp
     quadrature/GaussQuadrature.cpp
     quadrature/QuadratureBase.cpp
     quadrature/QuadratureRuleManager.cpp
     FiniteElementDiscretizationManager.cpp
     basis/BasisFunctionManager.cpp
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
                 HEADERS               ${finiteElement_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX finiteElement )

add_subdirectory( unitTests )

