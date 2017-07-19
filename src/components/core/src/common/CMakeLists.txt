
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	common/DataTypes.hpp
  	common/Logger.hpp
  	common/InterObjectRelation.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    common/DataTypes.cpp
    PARENT_SCOPE )
