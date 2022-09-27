# generate split test header file only if it does not yet exist
# in order to prevent a rebuild everytime cmake is configured
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h "")
  foreach(i RANGE 1 999)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/split_test_helper.h
      "#ifdef EIGEN_TEST_PART_${i}\n"
      "#define CALL_SUBTEST_${i}(FUNC) CALL_SUBTEST(FUNC)\n"
      "#else\n"
      "#define CALL_SUBTEST_${i}(FUNC)\n"
      "#endif\n\n"
    )
  endforeach()
endif()

set_property(GLOBAL PROPERTY EIGEN_CURRENT_SUBPROJECT "Unsupported")
