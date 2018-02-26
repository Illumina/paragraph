include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.7.0/include)

file(GLOB TEST_SOURCES test_*.cpp test_*.c)

add_executable(test_blackbox ${TEST_SOURCES})
target_link_libraries(test_blackbox ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS} ${GTEST_LIBRARY})
