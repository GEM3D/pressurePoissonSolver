add_library(MMWriter MMWriter.cpp)
target_link_libraries(MMWriter
    DomainCollection
    ${Trilinos_TPL_LIBRARIES}
    ${Trilinos_LIBRARIES} 
)

add_library(ClawWriter ClawWriter.cpp)
target_link_libraries(ClawWriter
    DomainCollection
    ${Trilinos_TPL_LIBRARIES}
    ${Trilinos_LIBRARIES} 
)

find_package(VTK NO_MODULE)

IF(VTK_FOUND)
    SET(CMAKE_CXX_FLAGS "-DHAVE_VTK ${CMAKE_CXX_FLAGS}")
    MESSAGE("VTK found, enabling vtk output")
ELSE()
    MESSAGE("VTK not found, disabling vtk output")
ENDIF()

if(VTK_FOUND)
    target_link_libraries(heat VtkWriter)
    add_library(VtkWriter VtkWriter.cpp)
    link_directories(${VTK_LIBRARY_DIRS})
    include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
    target_link_libraries(VtkWriter
        DomainCollection
        ${Trilinos_TPL_LIBRARIES}
        ${Trilinos_LIBRARIES} 
       ${MPI_CXX_LIBRARIES} 
       ${VTK_LIBRARIES} 
    )
endif()
