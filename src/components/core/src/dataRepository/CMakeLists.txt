
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/KeyIndexT.hpp
    dataRepository/ManagedGroup.hpp
    dataRepository/MappedVector.hpp    
    dataRepository/ViewWrapper.hpp
    dataRepository/ViewWrapperBase.hpp
    PARENT_SCOPE )

if( USE_ATK )
set(geosx_headers
  ${geosx_headers}
    dataRepository/SidreWrapper.hpp
    PARENT_SCOPE )
endif()

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/ManagedGroup.cpp
    dataRepository/ViewWrapper.cpp
    dataRepository/ViewWrapperBase.cpp
    PARENT_SCOPE )

if( USE_ATK )
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/SidreWrapper.cpp
    PARENT_SCOPE )
endif()