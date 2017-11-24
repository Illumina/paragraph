# Top-level Cmake file for Graph Tools
#
# Author: Peter Krusche <pkrusche@illumina.com>
#
cmake_minimum_required (VERSION 2.8.5)
project (GRM)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(cxx)
include(configureFiles)
include(FindZLIB)
include(FindBZip2)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scratch)

execute_process(
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/make_dependencies.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/scratch
    RESULT_VARIABLE DEPENDENCIES_BUILT)

if(NOT (${DEPENDENCIES_BUILT} STREQUAL "0"))
    message( FATAL_ERROR "Dependency build failed, CMake will exit." )
endif(NOT (${DEPENDENCIES_BUILT} STREQUAL "0"))

if (IS_DIRECTORY $ENV{PROTOBUF_INSTALL_PATH})
    message( "Using pre-built protobuf from $ENV{PROTOBUF_INSTALL_PATH}")
    set(CMAKE_PREFIX_PATH $ENV{PROTOBUF_INSTALL_PATH})
    set(Protobuf_DIR $ENV{PROTOBUF_INSTALL_PATH})
else()
    message( "Using included protobuf from ${CMAKE_BINARY_DIR}/scratch/protobuf-install")
    set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/scratch/protobuf-install)
    set(Protobuf_DIR ${CMAKE_BINARY_DIR}/scratch/protobuf-install)
endif ()

find_package(Protobuf REQUIRED)

include_directories(${PROTOBUF_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/src/protobuf)
include_directories(${CMAKE_BINARY_DIR}/scratch/spdlog/include)
include_directories(${GRM_SOURCE_DIR}/external/stream)
include_directories(${GRM_SOURCE_DIR}/external/src/c++/include)
include_directories(${GRM_SOURCE_DIR}/external/jsoncpp)
include_directories(${GRM_SOURCE_DIR}/external/klib)
include_directories(${GRM_SOURCE_DIR}/external/gssw)
include_directories(${GRM_SOURCE_DIR}/external/intervaltree)

include_directories(${GRM_SOURCE_DIR}/src/c++/include)

include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googletest/include)
include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googlemock/include)
include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googletest)
include_directories(${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googlemock)
add_library(gtest STATIC ${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googletest/src/gtest-all.cc)
add_library(gmock STATIC ${CMAKE_BINARY_DIR}/scratch/googletest-release-1.8.0/googlemock/src/gmock-all.cc)

set(GTEST_LIBRARY gtest gmock)


SET(DATA_SRC "${CMAKE_SOURCE_DIR}/share")
SET(DATA_DST "${CMAKE_BINARY_DIR}/share")
add_custom_target(copy_data ALL COMMENT "Copying ${DATA_SRC} -> ${DATA_DST}")

add_custom_command(TARGET copy_data
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${DATA_SRC} ${DATA_DST}
)

find_library(HTSLIB_LIBRARY NAMES libhts.a
             HINTS "${CMAKE_BINARY_DIR}/scratch/htslib" NO_DEFAULT_PATH)

find_package (Threads REQUIRED)

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)
find_package(Boost 1.5 COMPONENTS program_options filesystem system REQUIRED)

include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/scratch/htslib)
include_directories(${ZLIB_INCLUDE_DIR})
include_directories(${BZIP2_INCLUDE_DIR})
# boost sometimes generates warnings; we won't patch them so let's disable them using SYSTEM
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

# make libraries first
add_subdirectory (src/protobuf)
add_subdirectory (src/c++/lib)

add_subdirectory (external)

set(GRM_ALL_LIBS ${COMMON_LIBRARY}
                 ${HTSLIB_LIBRARY}
                 ${EXT_LIBRARY}
                 ${PROTO_LIBRARY}
                 ${GTEST_LIBRARY}
                 ${PROTOBUF_LIBRARIES}
                 ${CMAKE_THREAD_LIBS_INIT}
                 ${Boost_LIBRARIES}
                 ${BZIP2_LIBRARIES}
                 ${ZLIB_LIBRARIES}
        )

include(GetGitRevisionDescription)
git_describe(GRM_VERSION)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/grm/Version.hh.in"
               "${CMAKE_BINARY_DIR}/include/Version.hh")


add_subdirectory (src/c++/main)
add_subdirectory (src/c++/test)
add_subdirectory (src/c++/test-blackbox)
add_subdirectory (src/python)

# install
install (DIRECTORY "${CMAKE_BINARY_DIR}/bin/" DESTINATION bin FILES_MATCHING PATTERN "*" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE)
install (DIRECTORY "${CMAKE_BINARY_DIR}/lib/" DESTINATION lib)
install (DIRECTORY "${CMAKE_BINARY_DIR}/share/" DESTINATION share)
