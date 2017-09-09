
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )
    
add_subdirectory(codingUtilities)
add_subdirectory(common)
add_subdirectory(constitutive)
add_subdirectory(dataRepository)
add_subdirectory(fileIO)
add_subdirectory(finiteElement)
#add_subdirectory(legacy)
add_subdirectory(math)
add_subdirectory(schema)
add_subdirectory(managers)
add_subdirectory(MeshUtilities)
add_subdirectory(MPI_Communications)
add_subdirectory(PhysicsSolvers)
add_subdirectory(python)


# Python config
set( extraComponentsListLink "")
if( ENABLE_PYTHON )
  set( extraComponentsListLink ${extraComponentsListLink} python_interp )
endif()


if( ENABLE_MATHPRESSO )
  set( extraComponentsListLink ${extraComponentsListLink} mathpresso )
endif()


message("adding geosx_core library")
blt_add_library( NAME                  geosx_core
                 SOURCES               ${geosx_core_sources}
                 HEADERS               ${geosx_headers}
                 HEADERS_OUTPUT_SUBDIR geosx_core
                 DEPENDS_ON            cxx-utilities slic sidre conduit hdf5 pugixml spio silo trilinos optionparser ${extraComponentsListLink} ${thirdPartyLibs} )
                 
