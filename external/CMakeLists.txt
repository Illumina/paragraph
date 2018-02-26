
# External projects compiled into object library
set(LIBRARY_SOURCES gssw/gssw.c klib/ksw.c jsoncpp/jsoncpp.cpp)
set_source_files_properties(klib/ksw.c PROPERTIES COMPILE_FLAGS -msse4.1)
set_source_files_properties(gssw/gssw.c PROPERTIES COMPILE_FLAGS -msse4.1)
# no warnings. not ideal, but -Wno-enum-compare doesn't seem to do the trick
set_source_files_properties(jsoncpp/jsoncpp.cpp PROPERTIES COMPILE_FLAGS "-w")
add_library(external OBJECT ${LIBRARY_SOURCES})
set (EXT_LIBRARY external PARENT_SCOPE)


