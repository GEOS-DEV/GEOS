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
    PARENT_SCOPE )
