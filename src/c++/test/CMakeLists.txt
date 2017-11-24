
include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.7.0/include)

file(GLOB TEST_SOURCES test_*.cpp test_*.c)

add_executable(test_grm ${TEST_SOURCES})
target_link_libraries(test_grm ${GTEST_LIBRARY} ${GRM_LIBRARY} ${GRM_ALL_LIBS})
SET(LINK_LIBRARIES test_grm)
