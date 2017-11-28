project(core)

message( "Entering /src/components/core/CMakeLists.txt")

add_subdirectory(src)
add_subdirectory(tests)

#blt_add_code_check_targets(uncrustify.cfg)
geosx_add_code_checks(PREFIX core)

#install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)

message( "Leaving /src/components/core/CMakeLists.txt")
