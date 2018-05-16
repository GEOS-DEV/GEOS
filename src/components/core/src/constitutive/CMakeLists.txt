#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     LinearElasticIsotropic.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     LinearElasticIsotropic.cpp
   )

#set( dependencyList dataRepository )
#set( dependencyList2 "" )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  constitutive
                 SOURCES               ${constitutive_sources}
                 HEADERS               ${constitutive_headers}
                 DEPENDS_ON            dataRepository
                 SHARED FALSE
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX constitutive )