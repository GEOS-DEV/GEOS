configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_WaveEquationAcousticSolverAdjoint.xml
  ${CMAKE_CURRENT_BINARY_DIR}/test_WaveEquationAcousticSolverAdjoint.xml COPYONLY)
# This macro parse a test file and extract test_ function from it to pass it to ctest
macro(parse_test_file filename)
  set(filename ${filename})
  file(READ "${filename}" contents)
  string(REPLACE "\n" ";" rows ${contents})
  foreach(line ${rows})
    string(REGEX MATCHALL "^ *def *test_([A-Za-z_0-9])+ *" found ${line})
    if (found)
      string(REGEX REPLACE "^ *def *test_(([A-Za-z_0-9])+) *" "\\1" name "${found}")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pytest_template.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/${name}.sh
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
      add_test(NAME ${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${name}.sh)
    endif()
  endforeach()
endmacro(parse_test_file)

parse_test_file(test_WaveEquationAdjoint.py)
