
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/intrinsic/Wrapper.hpp
    dataRepository/intrinsic/WrapperBase.hpp
    dataRepository/intrinsic/WrapperCollection.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/intrinsic/Wrapper.cpp
    dataRepository/intrinsic/WrapperBase.cpp
    dataRepository/intrinsic/WrapperCollection.cpp
    PARENT_SCOPE )
