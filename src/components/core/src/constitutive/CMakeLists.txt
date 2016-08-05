#
# Specify all headers
#
set(geosx_core_sources
  ${geosx_core_sources}
    constitutive/ConstitutiveManager.hpp
    constitutive/ConstitutiveBase.hpp
    PARENT_SCOPE )
#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    constitutive/ConstitutiveManager.cpp
    constitutive/ConstitutiveBase.cpp
    PARENT_SCOPE )
