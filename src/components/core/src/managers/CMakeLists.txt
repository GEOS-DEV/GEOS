
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    managers/CommandLineParser.hpp
    managers/DomainPartition.hpp
    managers/EventManager.hpp
    managers/Functions/FunctionBase.hpp
    managers/Functions/SymbolicFunction.hpp
    managers/Functions/TableFunction.hpp
    managers/Functions/CompositeFunction.hpp
    managers/Functions/NewFunctionManager.hpp
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
    managers/EventManager.cpp
    managers/Functions/FunctionBase.cpp
    managers/Functions/SymbolicFunction.cpp
    managers/Functions/TableFunction.cpp
    managers/Functions/CompositeFunction.cpp
    managers/Functions/NewFunctionManager.cpp
    managers/ObjectManagerBase.cpp
    managers/ProblemManager.cpp    
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
    managers/ElementManager.cpp
    managers/ElementRegionT.cpp
    managers/FaceManager.cpp
    PARENT_SCOPE )
endif()