
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/CellBlockManager.hpp
    managers/CellBlock.hpp
    managers/CellBlockSubRegion.hpp
    managers/Event.hpp
    managers/EventManager.hpp
    managers/FaceManager.hpp
    managers/NodeManager.hpp
    managers/Functions/FunctionBase.hpp
    managers/Functions/SymbolicFunction.hpp
    managers/Functions/TableFunction.hpp
    managers/Functions/NewFunctionManager.hpp
    managers/ObjectManagerBase.hpp
    managers/ProblemManager.hpp
    managers/TableManager.hpp
    managers/ElementRegion.hpp
    managers/ElementRegionManager.hpp
    PARENT_SCOPE )


#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    managers/CommandLineParser.cpp
  	managers/DomainPartition.cpp
    managers/CellBlockManager.cpp
    managers/CellBlock.cpp
    managers/CellBlockSubRegion.cpp
    managers/EventManager.cpp
    managers/FaceManager.cpp
    managers/NodeManager.cpp
    managers/Functions/FunctionBase.cpp
    managers/Functions/SymbolicFunction.cpp
    managers/Functions/TableFunction.cpp
    managers/Functions/NewFunctionManager.cpp
    managers/ObjectManagerBase.cpp
    managers/ProblemManager.cpp
    managers/TableManager.cpp
    managers/ElementRegion.cpp
    managers/ElementRegionManager.cpp
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