cmake_minimum_required(VERSION 3.10)

add_executable(testPolyhedralDiscretizations testPolyhedralDiscretizations.cpp)

# Link with required libraries (replace coreComponents with actual target if needed)
target_link_libraries(testPolyhedralDiscretizations PRIVATE coreComponents)

add_test(NAME testPolyhedralDiscretizations COMMAND testPolyhedralDiscretizations)
