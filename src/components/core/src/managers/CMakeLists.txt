
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/ProblemManager.hpp
    managers/FaceManagerT.h
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/CommandLineParser.cpp
  	managers/DomainPartition.cpp
    managers/ProblemManager.cpp
    managers/FaceManagerT.cpp
    PARENT_SCOPE )
