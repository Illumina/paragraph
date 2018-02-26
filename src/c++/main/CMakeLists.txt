add_executable(grmpy grmpy.cpp)
target_link_libraries(grmpy ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})

add_executable(paragraph paragraph.cpp)
target_link_libraries(paragraph ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})

add_executable(kmerstats kmerstats.cpp)
target_link_libraries(kmerstats ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})

add_executable(idxdepth idxdepth.cpp)
target_link_libraries(idxdepth ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})

add_executable(pam pam.cpp)
target_link_libraries(pam ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})

add_executable(graph-to-fasta graph-to-fasta.cpp)
target_link_libraries(graph-to-fasta ${GRM_LIBRARY} ${GRM_EXTERNAL_LIBS})
