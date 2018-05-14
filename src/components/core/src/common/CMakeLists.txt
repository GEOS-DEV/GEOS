
#
# Specify all headers
#
set(common_headers
    integer_conversion.hpp
    DataTypes.hpp
  	InterObjectRelation.hpp
    Logger.hpp
    SortedArray.hpp )

#
# Specify all sources
#
set(common_sources
    Logger.cpp )

    
    
blt_add_library( NAME                  common
                 SOURCES               ${common_sources}
                 HEADERS               ${common_headers}
                 DEPENDS_ON            cxx-utilities mpi
                 SHARED FALSE
               )
               
               
geosx_add_code_checks(PREFIX common )