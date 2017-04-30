#
# Specify all headers
#
set(geosx_core_sources
  ${geosx_core_sources}
    finiteElement/FiniteElementSpace.hpp
    finiteElement/FiniteElementManager.hpp
    finiteElement/basis/BasisBase.hpp
    finiteElement/basis/LagrangeBasis.hpp
    finiteElement/basis/Polynomial.hpp    
    finiteElement/quadrature/QuadratureBase.hpp
    finiteElement/quadrature/GaussQuadrature.hpp
    PARENT_SCOPE )
#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    finiteElement/FiniteElementSpace.cpp
    finiteElement/FiniteElementManager.cpp
    finiteElement/basis/BasisBase.cpp
    finiteElement/basis/LagrangeBasis.cpp
    finiteElement/basis/Polynomial.cpp    
    finiteElement/quadrature/QuadratureBase.cpp
    finiteElement/quadrature/GaussQuadrature.cpp
    PARENT_SCOPE )
