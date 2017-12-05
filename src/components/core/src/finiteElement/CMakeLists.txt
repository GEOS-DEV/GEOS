#
# Specify all headers
#
set(geosx_core_sources
  ${geosx_core_sources}
    finiteElement/FiniteElementSpace.hpp
    finiteElement/FiniteElementManager.hpp
    finiteElement/ElementLibrary/FiniteElement.h
    finiteElement/ElementLibrary/FiniteElementBase.h
    finiteElement/Kinematics.h
    finiteElement/basis/BasisBase.hpp
    finiteElement/basis/LagrangeBasis.hpp
    finiteElement/basis/Polynomial.hpp    
    finiteElement/quadrature/QuadratureBase.hpp
    finiteElement/quadrature/GaussQuadrature.hpp
    finiteElement/FiniteElementSpaceManager.hpp
    finiteElement/basis/BasisFunctionManager.hpp
    finiteElement/quadrature/QuadratureRuleManager.hpp
    PARENT_SCOPE )
#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    finiteElement/FiniteElementSpace.cpp
    finiteElement/FiniteElementManager.cpp
    finiteElement/ElementLibrary/FiniteElement.cpp
    finiteElement/ElementLibrary/FiniteElementBase.cpp
    finiteElement/Kinematics.cpp
    finiteElement/basis/BasisBase.cpp
    finiteElement/basis/LagrangeBasis.cpp
    finiteElement/basis/Polynomial.cpp    
    finiteElement/quadrature/QuadratureBase.cpp
    finiteElement/quadrature/GaussQuadrature.cpp
    finiteElement/FiniteElementSpaceManager.cpp
    finiteElement/basis/BasisFunctionManager.cpp
    finiteElement/quadrature/QuadratureRuleManager.cpp
    PARENT_SCOPE )
