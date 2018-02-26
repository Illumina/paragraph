##
## file CMakeLists.txt
##
## Configuration file for the validation subfolder
##
## author Roman Petrovski
##

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/simulate-reads.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/simulate-reads.sh" @ONLY)
install (FILES "${CMAKE_CURRENT_BINARY_DIR}/simulate-reads.sh" DESTINATION libexec PERMISSIONS ${PARAGRAPH_EXEC_PERMISSIONS})
install (DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION libexec FILES_MATCHING PATTERN "*.sh" PERMISSIONS ${PARAGRAPH_EXEC_PERMISSIONS})
