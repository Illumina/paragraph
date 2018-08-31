# Top-level Cmake file for Graph Tools
#
# Author: Peter Krusche <pkrusche@illumina.com>
#
cmake_minimum_required (VERSION 2.8.5)
project (GRM)

# Support for static linking.
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
if    (PARAGRAPH_FORCE_STATIC_LINK)
    message(STATUS "All libraries will be statically linked")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
endif (PARAGRAPH_FORCE_STATIC_LINK)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(cxx)
include(configureFiles)

include(FindZLIB)
include(GetBoost)
include(GetHtslib)
include(GetGoogleTest)
include(GetGraphTools)
include(GetSpdlog)

include_directories(${GRM_SOURCE_DIR}/external/src/c++/include)
include_directories(${GRM_SOURCE_DIR}/external/jsoncpp)
include_directories(${GRM_SOURCE_DIR}/external/klib)
include_directories(${GRM_SOURCE_DIR}/external/gssw)
include_directories(${GRM_SOURCE_DIR}/src/c++/include)

SET(DATA_SRC "${CMAKE_SOURCE_DIR}/share")
SET(DATA_DST "${CMAKE_BINARY_DIR}/share")
add_custom_target(copy_data ALL COMMENT "Copying ${DATA_SRC} -> ${DATA_DST}")

add_custom_command(TARGET copy_data
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${DATA_SRC} ${DATA_DST}
)

find_package (Threads REQUIRED)

include_directories(${ZLIB_INCLUDE_DIR})
include_directories(${BZIP2_INCLUDE_DIR})

# make libraries first
add_subdirectory (src/c++/lib)
add_subdirectory (external)

add_library(grm STATIC $<TARGET_OBJECTS:grmpy_common> $<TARGET_OBJECTS:external>)
set(GRM_LIBRARY grm)

set(GRM_EXTERNAL_LIBS ${HTSLIB_LIBRARY}
                      ${GRAPHTOOLS_LIBRARY}
                      ${CMAKE_THREAD_LIBS_INIT}
                      ${Boost_LIBRARIES}
                      ${BZIP2_LIBRARIES}
                      ${ZLIB_LIBRARIES}
        )

set(PARAGRAPH_EXEC_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE)

include(GetGitRevisionDescription)
git_describe(GRM_VERSION)
if ("${GRM_VERSION}" STREQUAL "HEAD-HASH-NOTFOUND")
SET(PARAGRAPH_VERSION_MAJOR "MAJOR")
SET(PARAGRAPH_VERSION_MINOR "MINOR")
SET(PARAGRAPH_VERSION_PATCH "PATCH")
else("${GRM_VERSION}" STREQUAL "HEAD-HASH-NOTFOUND")
STRING(REGEX REPLACE "^paragraph-v([0-9]+)\\.[0-9]+\\-.+" "\\1" PARAGRAPH_VERSION_MAJOR "${GRM_VERSION}")
STRING(REGEX REPLACE "^paragraph-v[0-9]+\\.([0-9])+\\-.+" "\\1" PARAGRAPH_VERSION_MINOR "${GRM_VERSION}")
STRING(REGEX REPLACE "^paragraph-v[0-9]+\\.[0-9]+\\-(.+)" "\\1" PARAGRAPH_VERSION_PATCH "${GRM_VERSION}")
endif ("${GRM_VERSION}" STREQUAL "HEAD-HASH-NOTFOUND")

if (CPACK_GENERATOR)
    message (STATUS "Configuring to produce the following package types: ${CPACK_GENERATOR}")
    SET(CPACK_PACKAGE_VENDOR "Illumina")
    SET(CPACK_PACKAGE_VERSION_MAJOR "${PARAGRAPH_VERSION_MAJOR}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${PARAGRAPH_VERSION_MINOR}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${PARAGRAPH_VERSION_PATCH}")
    INCLUDE(CPack)
endif (CPACK_GENERATOR)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/grm/Version.hh.in"
               "${CMAKE_BINARY_DIR}/include/grm/Version.hh")
include_directories("${CMAKE_BINARY_DIR}/include")

add_subdirectory (src/c++/main)
add_subdirectory (src/c++/test)
add_subdirectory (src/c++/test-blackbox)
add_subdirectory (src/python)
add_subdirectory (src/make/validation)
add_subdirectory (src/sh/validation)

# install binaries + libraries
install (DIRECTORY "${CMAKE_BINARY_DIR}/bin/" DESTINATION bin FILES_MATCHING PATTERN "*" PERMISSIONS ${PARAGRAPH_EXEC_PERMISSIONS})
install (TARGETS grm DESTINATION lib)

# install headers
install (DIRECTORY "${CMAKE_SOURCE_DIR}/src/c++/include/" DESTINATION include)
install (FILES "${CMAKE_BINARY_DIR}/include/grm/Version.hh" DESTINATION include/grm)
install (DIRECTORY "${CMAKE_SOURCE_DIR}/external/gssw" DESTINATION include FILES_MATCHING PATTERN "*.h")
install (DIRECTORY "${CMAKE_SOURCE_DIR}/external/jsoncpp/json" DESTINATION include)

# install docs
install (DIRECTORY "${CMAKE_SOURCE_DIR}/doc/" DESTINATION doc)
install (DIRECTORY "${CMAKE_BINARY_DIR}/share/" DESTINATION share)
install (DIRECTORY "${CMAKE_BINARY_DIR}/lib/python3" DESTINATION lib)
