add_library(MMWriter MMWriter.cpp)
target_link_libraries(MMWriter
    ${Trilinos_TPL_LIBRARIES}
    ${Trilinos_LIBRARIES} 
)

add_library(ClawWriter ClawWriter.cpp)
target_link_libraries(ClawWriter
    ${Trilinos_TPL_LIBRARIES}
    ${Trilinos_LIBRARIES} 
)


if(VTK_FOUND)
    add_library(VtkWriter VtkWriter.cpp VtkWriter2d.cpp)
    link_directories(${VTK_LIBRARY_DIRS})
    include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
    target_link_libraries(VtkWriter
        ${Trilinos_TPL_LIBRARIES}
        ${Trilinos_LIBRARIES} 
       ${MPI_CXX_LIBRARIES} 
       ${VTK_LIBRARIES} 
    )
endif()
