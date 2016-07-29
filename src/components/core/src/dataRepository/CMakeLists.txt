
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/SidreWrapper.hpp
    dataRepository/Wrapper.hpp
    dataRepository/WrapperBase.hpp
    dataRepository/WrapperCollection.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/SidreWrapper.cpp
    dataRepository/Wrapper.cpp
    dataRepository/WrapperBase.cpp
    dataRepository/WrapperCollection.cpp
    PARENT_SCOPE )
