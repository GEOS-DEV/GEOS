
#blt_add_component(COMPONENT_NAME optionparser DEFAULT_STATE ON)
add_subdirectory(optionparser)
add_subdirectory(pugixml)

# muparserx settings
set(ENABLE_MUPARSERX ON CACHE BOOL  "Enables muparserx Plugin")
if( ENABLE_MUPARSERX )
  message("adding muparserx library")
  set(BUILD_EXAMPLES OFF CACHE BOOL "Disables muparserx example")
  add_subdirectory(muparserx)
endif()
