#
# Specify all headers
#
set( coupling_headers
     TribolCoupling.hpp
   )
#
# Specify all sources
#
set( coupling_sources
     TribolCoupling.cpp
   )

blt_add_library( NAME                  coupling
                 SOURCES               ${coupling_sources}
                 HEADERS               ${coupling_headers}
                 DEPENDS_ON            dataRepository mpi tribol worlds_base
               )
               
target_include_directories( coupling PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )
geosx_add_code_checks( PREFIX coupling )
