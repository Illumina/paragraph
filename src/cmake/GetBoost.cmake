
if (IS_DIRECTORY ${BOOST_ROOT})
    message( "Using pre-built boost from ${BOOST_ROOT}")
else()
    message( "Building included subset of Boost" )
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
    include(ExternalProject) \n
    ExternalProject_Add(boost \n
      URL \"${CMAKE_SOURCE_DIR}/external/boost_subset_1_67_0.tar.gz\" \n
      URL_MD5 9b6dce185b01cd34c85ec020de805a9b \n
      BUILD_IN_SOURCE 1 \n
      UPDATE_COMMAND \"\" \n
      PATCH_COMMAND \"\" \n
      SOURCE_DIR \"${CMAKE_BINARY_DIR}/external/boost-src\" \n
      CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_COMMAND} \n
      BUILD_COMMAND  ${BOOST_B2_COMMAND} install \n
        --prefix=${BOOST_ROOT} \n
        --threading=single,multi \n
        --link=static \n
        --variant=${BOOST_VARIANT} \n
        -j4 \n
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
endif()

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)
find_package(Boost 1.58 COMPONENTS iostreams program_options filesystem system REQUIRED)

# boost sometimes generates warnings; we won't patch them so let's disable them using SYSTEM
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

