
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  fileIO/xmlWrapper.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
  fileIO/xmlWrapper.cpp
    PARENT_SCOPE )
