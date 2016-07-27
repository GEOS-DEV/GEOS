
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	fileIO/ticpp/HierarchicalDataNode.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    fileIO/ticpp/HierarchicalDataNode.cpp
    PARENT_SCOPE )

message( "CMAKE_CURRENT_LIST_DIR =${CMAKE_CURRENT_LIST_DIR}" )