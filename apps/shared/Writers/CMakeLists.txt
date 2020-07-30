#determine sources
list(APPEND UTILS_WRITERS_SRCS ClawWriter.cpp)
if(VTK_FOUND)
    list(APPEND UTILS_WRITERS_SRCS VtkWriter.cpp VtkWriter2d.cpp)
ENDIF()

add_library(UTILS_WRITERS ${UTILS_WRITERS_SRCS})

if(VTK_FOUND)
    link_directories(${VTK_LIBRARY_DIRS})
    target_include_directories(UTILS_WRITERS PUBLIC ${VTK_INCLUDE_DIRS})
endif()

target_link_libraries(UTILS_WRITERS Thunderegg)
if(VTK_FOUND)
    target_link_libraries(UTILS_WRITERS ${VTK_LIBRARIES})
endif()
