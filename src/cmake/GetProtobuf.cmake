set(PROTOBUF_INSTALL_PATH "$ENV{PROTOBUF_INSTALL_PATH}" CACHE STRING "Specify a pre-built version of protobuf (install prefix).")

# find or install protobuf
if (IS_DIRECTORY ${PROTOBUF_INSTALL_PATH})
    message( "Using pre-built protobuf from ${PROTOBUF_INSTALL_PATH}")
    set(Protobuf_DIR ${PROTOBUF_INSTALL_PATH})
else()
    message( "Using included protobuf from ${CMAKE_BINARY_DIR}/scratch/protobuf-install")
    set(Protobuf_DIR ${CMAKE_BINARY_DIR}/external/protobuf-install)
    set(PROTOBUF_INSTALL_PATH ${CMAKE_BINARY_DIR}/external/protobuf-install)

    FILE(WRITE "${CMAKE_BINARY_DIR}/external/protobuf-build/CMakeLists.txt" "\
        cmake_minimum_required(VERSION 2.8.5)
        project(protobuf-build NONE)
        include(ExternalProject)
        ExternalProject_Add(protobuf
        URL               \"${CMAKE_SOURCE_DIR}/external/protobuf-cpp-3.3.0.tar.gz\"
        URL_HASH          MD5=73c28d3044e89782bdc8d9fdcfbb5792
        SOURCE_DIR        \"${CMAKE_BINARY_DIR}/external/protobuf-src\"
        INSTALL_DIR       \"${CMAKE_BINARY_DIR}/external/protobuf-install\"
        CONFIGURE_COMMAND cd <SOURCE_DIR> && ./configure --enable-shared=no --prefix=<INSTALL_DIR>
        BUILD_COMMAND make -C <SOURCE_DIR> -j8
        INSTALL_COMMAND make -C <SOURCE_DIR> install)")

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/protobuf-build )
    if(result)
        message(FATAL_ERROR "CMake step for protobuf failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/protobuf-build )
    if(result)
        message(FATAL_ERROR "Build step for protobuf failed: ${result}")
    endif()

endif ()

set(CMAKE_PREFIX_PATH ${PROTOBUF_INSTALL_PATH})
# set(Protobuf_DEBUG 1)
find_package(Protobuf 3.3.0 REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS})
