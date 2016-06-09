
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    codingUtilities/stackTrace.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    codingUtilities/stackTrace.cpp
    PARENT_SCOPE )
