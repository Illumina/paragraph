##
## file CMakeLists.txt
##
## Configuration file for the validation subfolder
##
## author Roman Petrovski
##

install (DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION libexec FILES_MATCHING PATTERN "*.mk")
