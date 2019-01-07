#
# Specify all headers
#
set( finiteElement_headers
     FiniteElementDiscretization.hpp
     ElementLibrary/FiniteElement.h
     ElementLibrary/FiniteElementBase.h
     Kinematics.h
     basis/BasisBase.hpp
     basis/LagrangeBasis.hpp
     basis/Polynomial.hpp    
     quadrature/QuadratureBase.hpp
     quadrature/GaussQuadrature.hpp
     FiniteElementDiscretizationManager.hpp
     basis/BasisFunctionManager.hpp
     quadrature/QuadratureRuleManager.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     ElementLibrary/FiniteElement.cpp
     ElementLibrary/FiniteElementBase.cpp
     Kinematics.cpp
     basis/BasisBase.cpp
     basis/LagrangeBasis.cpp
     basis/Polynomial.cpp    
     quadrature/QuadratureBase.cpp
     quadrature/GaussQuadrature.cpp
     FiniteElementDiscretizationManager.cpp
     basis/BasisFunctionManager.cpp
     quadrature/QuadratureRuleManager.cpp
   )

blt_add_library( NAME                  finiteElement
                 SOURCES               ${finiteElement_sources}
                 HEADERS               ${finiteElement_headers}
                 DEPENDS_ON            dataRepository codingUtilities
               )
               
target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX finiteElement )