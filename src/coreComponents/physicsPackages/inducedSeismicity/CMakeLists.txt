# Specify solver headers
set( physicsPackages_headers
     ${physicsPackages_headers}
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SeismicityRateKernels.hpp 
     PARENT_SCOPE )

# Specify solver sources
set( physicsPackages_sources
     ${physicsPackages_sources}
     inducedSeismicity/SeismicityRate.cpp 
     PARENT_SCOPE )
