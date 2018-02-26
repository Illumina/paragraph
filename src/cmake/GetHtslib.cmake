set(HTSLIB_INSTALL_PATH "$ENV{HTSLIB_INSTALL_PATH}" CACHE STRING "Specify a pre-built version of htslib (install prefix).")

include(FindZLIB)
include(FindBZip2)
include(FindLibLZMA)
if (IS_DIRECTORY ${HTSLIB_INSTALL_PATH})
    message( "Using pre-built htslib from ${HTSLIB_INSTALL_PATH}")
else()
    message( "Using included htslib" )
    set(HTSLIB_INSTALL_PATH ${CMAKE_BINARY_DIR}/external/htslib-install)

    FILE(WRITE "${CMAKE_BINARY_DIR}/external/htslib-build/CMakeLists.txt" "
        cmake_minimum_required(VERSION 2.8.5)\n
        project(htslib-build NONE)
        include(ExternalProject)
        ExternalProject_Add(htslib
        URL               \"${CMAKE_SOURCE_DIR}/external/htslib.tar.gz\"
        URL_HASH          MD5=fde78c5212a41921f68665617b03f371
        SOURCE_DIR        \"${CMAKE_BINARY_DIR}/external/htslib-src\"
        INSTALL_DIR       \"${HTSLIB_INSTALL_PATH}\"
        CONFIGURE_COMMAND \"\"
        BUILD_COMMAND make -C <SOURCE_DIR> prefix=<INSTALL_DIR>
        INSTALL_COMMAND make -C <SOURCE_DIR> install prefix=<INSTALL_DIR> )")

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/htslib-build )
    if(result)
        message(FATAL_ERROR "CMake step for htslib failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/htslib-build )
    if(result)
        message(FATAL_ERROR "Build step for htslib failed: ${result}")
    endif()
endif()

# find library + add include path
include_directories(${ZLIB_INCLUDE_DIR})
include_directories(${BZIP2_INCLUDE_DIR})
include_directories(${LIBLZMA_INCLUDE_DIRS})

include_directories("${HTSLIB_INSTALL_PATH}/include")
find_library(HTSLIB_LIBRARY NAMES libhts.a
        HINTS "${HTSLIB_INSTALL_PATH}/lib" NO_DEFAULT_PATH)
set(HTSLIB_LIBRARY ${HTSLIB_LIBRARY} ${LIBLZMA_LIBRARIES})

# install htslib
FILE(GLOB HTSLIB_LIBRARY_FILES "${HTSLIB_INSTALL_PATH}/lib/*.a" "${HTSLIB_INSTALL_PATH}/lib/*.dylib" "${HTSLIB_INSTALL_PATH}/lib/*.so")
install (FILES ${HTSLIB_LIBRARY_FILES} DESTINATION lib)
install (DIRECTORY "${HTSLIB_INSTALL_PATH}/include/htslib" DESTINATION "include")
