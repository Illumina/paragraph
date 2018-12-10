
# BOOST_ROOT > System > included
set(USE_SYSTEM_BOOST TRUE CACHE BOOL "Enable/disable the use of system-wide Boost installation")

if (IS_DIRECTORY ${BOOST_ROOT})
    message( "Using pre-built boost from ${BOOST_ROOT}")
elseif (${USE_SYSTEM_BOOST})
    find_package(Boost 1.58 COMPONENTS iostreams program_options filesystem system REQUIRED)
else()
    message( "Building included Boost" )
    set( BOOST_BOOTSTRAP_COMMAND )
    if( UNIX )
      set( BOOST_BOOTSTRAP_COMMAND ./bootstrap.sh )
      set( BOOST_B2_COMMAND ./b2 )
    else()
      if( WIN32 )
        set( BOOST_BOOTSTRAP_COMMAND bootstrap.bat )
        set( BOOST_B2_COMMAND b2.exe )
      endif()
    endif()

    set( BOOST_ROOT ${CMAKE_BINARY_DIR}/external/boost-install )

    set (BOOST_VARIANT "release")
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set (BOOST_VARIANT "debug")
    endif()

    FILE(WRITE "${CMAKE_BINARY_DIR}/external/boost-build/CMakeLists.txt" "
    cmake_minimum_required(VERSION 3.1.0)
    include(ExternalProject) \n
    ExternalProject_Add(boost \n
      URL \"https://downloads.sourceforge.net/project/boost/boost/1.61.0/boost_1_61_0.tar.gz\" \n
      URL_MD5 874805ba2e2ee415b1877ef3297bf8ad \n
      BUILD_IN_SOURCE 1 \n
      UPDATE_COMMAND \"\" \n
      PATCH_COMMAND \"\" \n
      SOURCE_DIR \"${CMAKE_BINARY_DIR}/external/boost-src\" \n
      CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_COMMAND} \n
      BUILD_COMMAND  ${BOOST_B2_COMMAND} install \n
        --prefix=${BOOST_ROOT} \n
        --threading=multi \n
        --link=static \n
        --variant=${BOOST_VARIANT} \n
        -j4 > ${CMAKE_BINARY_DIR}/boost_build.log \n
      INSTALL_COMMAND \"\" \n
      INSTALL_DIR \"\"\n
    )")

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/boost-build )
    if(result)
        message(FATAL_ERROR "CMake step for boost failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/boost-build )
    if(result)
        message(FATAL_ERROR "Build step for boost failed: ${result}")
    endif()

    set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
    set(Boost_USE_MULTITHREADED      ON)
    set(Boost_USE_STATIC_RUNTIME     ON)
endif()

find_package(Boost 1.58 COMPONENTS iostreams program_options filesystem system REQUIRED)

# boost sometimes generates warnings; we won't patch them so let's disable them using SYSTEM
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
