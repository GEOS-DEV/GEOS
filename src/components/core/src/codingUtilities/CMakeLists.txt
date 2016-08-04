
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    codingUtilities/StringUtilities.hpp
    codingUtilities/UnitManager.h
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    codingUtilities/StringUtilities.cpp
    codingUtilities/UnitManager.cpp
    PARENT_SCOPE )
