
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/SidreWrapper.hpp
    dataRepository/SynchronizedGroup.hpp
    dataRepository/WrapperView.hpp
    dataRepository/WrapperViewBase.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/SidreWrapper.cpp
    dataRepository/SynchronizedGroup.cpp
    dataRepository/WrapperView.cpp
    dataRepository/WrapperViewBase.cpp
    PARENT_SCOPE )
