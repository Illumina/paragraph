file(GLOB_RECURSE GRM_LIBRARY_SOURCES Graph.cpp graphs/*.cpp genotyping/*.cpp paragraph/*.cpp grm/*.cpp grmpy/*.cpp
        common/*.cpp variant/*.cpp intervals/*.cpp idxdepth/*.cpp pam/*.cpp)
add_library(grmpy_common STATIC ${GRM_LIBRARY_SOURCES})
set(GRM_LIBRARY grmpy_common PARENT_SCOPE)

