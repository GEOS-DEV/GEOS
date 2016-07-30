
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	managers/ProblemManager.hpp
  	managers/DomainPartition.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/ProblemManager.cpp
  	managers/DomainPartition.cpp
    PARENT_SCOPE )
