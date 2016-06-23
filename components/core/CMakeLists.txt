project(core)

add_subdirectory(src)
add_subdirectory(tests)

add_code_check_targets(uncrustify.cfg)

#install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
