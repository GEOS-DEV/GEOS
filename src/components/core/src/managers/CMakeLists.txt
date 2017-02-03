
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/ElementManager.hpp
    managers/ElementRegion.hpp
    managers/EventManager.hpp
    managers/FaceManager.hpp
    managers/NodeManager.hpp
    managers/ObjectManagerBase.hpp
    managers/ProblemManager.hpp
    managers/TableManager.hpp
    PARENT_SCOPE )


#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/CommandLineParser.cpp
  	managers/DomainPartition.cpp
    managers/ElementManager.cpp
    managers/ElementRegion.cpp
    managers/EventManager.cpp
    managers/FaceManager.cpp
    managers/NodeManager.cpp
    managers/ObjectManagerBase.cpp
    managers/ProblemManager.cpp
    managers/TableManager.cpp
    PARENT_SCOPE )

if( 0 )
set(geosx_headers
  ${geosx_headers}
    managers/EdgeManager.hpp
    managers/ElementManager.hpp
    managers/ElementRegionT.hpp
    managers/FaceManager.hpp
    PARENT_SCOPE )
    
set(geosx_core_sources
  ${geosx_core_sources}
    managers/EdgeManager.cpp
    managers/ElementManager.cpp
    managers/ElementRegionT.cpp
    managers/FaceManager.cpp
    PARENT_SCOPE )
endif()