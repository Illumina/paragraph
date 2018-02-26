######################### Google Test ############################
# Download and unpack spdlog at configure time

FILE(WRITE "${CMAKE_BINARY_DIR}/external/spdlog-build/CMakeLists.txt" "\
cmake_minimum_required(VERSION 2.8.5)
project(spdlog-build NONE)
include(ExternalProject)
ExternalProject_Add(spdlog
        URL               \"${CMAKE_SOURCE_DIR}/external/spdlog.tar.gz\"
        URL_HASH          MD5=8c0faa7d7aa02500b2a3e61563bfed75
        SOURCE_DIR        \"${CMAKE_BINARY_DIR}/external/spdlog-src\"
        BINARY_DIR        \"${CMAKE_BINARY_DIR}/external/spdlog-build\"
        CONFIGURE_COMMAND \"\"
        BUILD_COMMAND     \"\"
        INSTALL_COMMAND   \"\"
        TEST_COMMAND      \"\"
        )"
        )

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/spdlog-build )
if(result)
    message(FATAL_ERROR "CMake step for spdlog failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/spdlog-build )
if(result)
    message(FATAL_ERROR "Build step for spdlog failed: ${result}")
endif()

include_directories("${CMAKE_BINARY_DIR}/external/spdlog-src/include")
# bundle headers for installation
install (DIRECTORY "${CMAKE_BINARY_DIR}/external/spdlog-src/include/" DESTINATION include)
##################################################################
