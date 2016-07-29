
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	codingUtilities/macros.hpp
    codingUtilities/ObjectCatalog.hpp
    codingUtilities/stackTrace.hpp
    codingUtilities/SetSignalHandling.hpp
    codingUtilities/sfinae.hpp
    codingUtilities/StringUtilities.hpp
    codingUtilities/UnitManager.h
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    codingUtilities/stackTrace.cpp
    codingUtilities/SetSignalHandling.cpp
    codingUtilities/StringUtilities.cpp
    codingUtilities/UnitManager.cpp
    PARENT_SCOPE )
