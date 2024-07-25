# Specify solver headers
set( physicsPackages_headers
      ${physicsPackages_headers}
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp 
     PARENT_SCOPE )

# Specify solver sources
set( physicsPackages_sources
      ${physicsPackages_sources}
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp 
     PARENT_SCOPE )
