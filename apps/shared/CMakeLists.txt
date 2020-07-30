#determine sources first
list(APPEND UTILS_SRCS Init.cpp Init.h)
if(p4est_FOUND)
    list(APPEND UTILS_SRCS TreeToP4est.cpp TreeToP4est.h)
endif()

add_library(UTILS ${UTILS_SRCS})

target_include_directories(UTILS PUBLIC ".")

target_link_libraries(UTILS Thunderegg)
target_link_libraries(UTILS UTILS_WRITERS)

add_subdirectory(Writers)
