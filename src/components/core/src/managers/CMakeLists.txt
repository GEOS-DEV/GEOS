
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/ProblemManager.hpp
    managers/ElementManagerT.hpp
    managers/FaceManager.hpp
    managers/ObjectManagerBase.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/CommandLineParser.cpp
  	managers/DomainPartition.cpp
    managers/ProblemManager.cpp
    managers/ElementManagerT.cpp
    managers/FaceManager.cpp
    managers/ObjectManagerBase.cpp
    PARENT_SCOPE )
