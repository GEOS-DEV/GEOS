
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	array/MultidimensionalArray.hpp
    PARENT_SCOPE )

#
# Specify all sources
#

blt_add_executable(NAME array
                   SOURCES main.cpp
                   DEPENDS_ON ${GEOSX_LINK_PREPEND_FLAG} geosx_core ${GEOSX_LINK_POSTPEND_FLAG}
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja 
                   )
