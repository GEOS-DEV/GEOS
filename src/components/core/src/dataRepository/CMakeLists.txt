
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/KeyIndexT.hpp
    dataRepository/ManagedGroup.hpp
    dataRepository/MappedVector.hpp    
    dataRepository/SidreWrapper.hpp
    dataRepository/ViewWrapper.hpp
    dataRepository/ViewWrapperBase.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/ManagedGroup.cpp
    dataRepository/SidreWrapper.cpp
    dataRepository/ViewWrapper.cpp
    dataRepository/ViewWrapperBase.cpp
    PARENT_SCOPE )

