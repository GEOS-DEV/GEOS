
#blt_add_component(COMPONENT_NAME optionparser DEFAULT_STATE ON)
add_subdirectory(optionparser)
add_subdirectory(pugixml)

# muparserx settings
if( ENABLE_MUPARSERX )
  message("adding muparserx library")
  set(BUILD_EXAMPLES OFF CACHE BOOL "Disables muparserx example")
  add_subdirectory(muparserx)
endif()
