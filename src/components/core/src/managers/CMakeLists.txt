
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/EdgeManager.hpp
    managers/ElementManager.hpp
    managers/FaceManager.hpp
    managers/ObjectManagerBase.hpp
    managers/ProblemManager.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/CommandLineParser.cpp
  	managers/DomainPartition.cpp
    managers/EdgeManager.cpp
    managers/ElementManager.cpp
    managers/FaceManager.cpp
    managers/ObjectManagerBase.cpp
    managers/ProblemManager.cpp
    PARENT_SCOPE )
