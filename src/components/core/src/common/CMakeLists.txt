
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    common/integer_conversion.hpp
    common/DataTypes.hpp
  	common/InterObjectRelation.hpp
    common/Logger.hpp
    common/SortedArray.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    common/Logger.cpp
    PARENT_SCOPE )
