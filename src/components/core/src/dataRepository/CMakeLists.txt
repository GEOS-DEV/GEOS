
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/SidreWrapper.hpp
    dataRepository/ManagedGroup.hpp
    dataRepository/ViewWrapper.hpp
    dataRepository/ViewWrapperBase.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/SidreWrapper.cpp
    dataRepository/ManagedGroup.cpp
    dataRepository/ViewWrapper.cpp
    dataRepository/ViewWrapperBase.cpp
    PARENT_SCOPE )
