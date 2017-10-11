
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  fileIO/xmlWrapper.hpp
  fileIO/silo/SiloFile.hpp
  fileIO/blueprint/Blueprint.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
  fileIO/xmlWrapper.cpp
  fileIO/silo/SiloFile.cpp
  fileIO/blueprint/Blueprint.cpp
    PARENT_SCOPE )
