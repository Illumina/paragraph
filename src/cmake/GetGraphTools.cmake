######################### Google Test ############################
# Download and unpack googletest at configure time


if(NOT DEFINED GRAPHTOOLS_SOURCE_DIR)
    set(GRAPHTOOLS_URL         "${CMAKE_SOURCE_DIR}/external/graph-tools.tar.gz")
    set(GRAPHTOOLS_URL_HASH    "MD5=b2e90e588a68b572a8b7bccc45e2e451")
    set(GRAPHTOOLS_SOURCE_DIR  "${CMAKE_BINARY_DIR}/external/graphtools-src")
endif (NOT DEFINED GRAPHTOOLS_SOURCE_DIR)

FILE(WRITE "${CMAKE_BINARY_DIR}/external/graphtools-build/CMakeLists.txt" "
cmake_minimum_required(VERSION 3.1.0)
project(graphtools-build NONE)
include(ExternalProject)

ExternalProject_Add(graphtools
        URL               \"${GRAPHTOOLS_URL}\"
        URL_HASH          \"${GRAPHTOOLS_URL_HASH}\"
        SOURCE_DIR        \"${GRAPHTOOLS_SOURCE_DIR}\"
        BINARY_DIR        \"${CMAKE_BINARY_DIR}/external/graphtools-build\"
        CONFIGURE_COMMAND \"\"
        BUILD_COMMAND     \"\"
        INSTALL_COMMAND   \"\"
        TEST_COMMAND      \"\"
        )"
)

function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

set(BOOST_OPTION "")

getListOfVarsStartingWith("BOOST_" matchedVars)
foreach (_var IN LISTS matchedVars)
    #message("Passing Boost option: ${_var}=${${_var}}")
    set(BOOST_OPTION " -D${_var}=${${_var}}")
endforeach()
getListOfVarsStartingWith("Boost_" matchedVars)
foreach (_var IN LISTS matchedVars)
    #message("Passing Boost option: ${_var}=${${_var}}")
    set(BOOST_OPTION " -D${_var}=${${_var}}")
endforeach()

execute_process(COMMAND ${CMAKE_COMMAND} ${BOOST_OPTION} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/graphtools-build )
if(result)
    message(FATAL_ERROR "CMake step for graphtools failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/graphtools-build )
if(result)
    message(FATAL_ERROR "Build step for graphtools failed: ${result}")
endif()

# Add graphtools directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${GRAPHTOOLS_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}/external/graphtools-build)

include_directories(${GRAPHTOOLS_SOURCE_DIR}/include)
SET(GRAPHTOOLS_LIBRARY graphtools)

install(DIRECTORY ${GRAPHTOOLS_SOURCE_DIR}/include/graphcore DESTINATION include/)
install(DIRECTORY ${GRAPHTOOLS_SOURCE_DIR}/include/graphalign DESTINATION include/)
install(DIRECTORY ${GRAPHTOOLS_SOURCE_DIR}/include/graphutils DESTINATION include/)
install(FILES $<TARGET_FILE:graphtools> DESTINATION lib/)

##################################################################
