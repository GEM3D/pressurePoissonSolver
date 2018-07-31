add_executable(test SchurDomain.cpp Domain.cpp GMG.cpp test.cpp Side.cpp Octant.cpp OctTree.cpp
    DomainCollection.cpp)
target_link_libraries(test
    ${MPI_CXX_LIBRARIES} 
    ${PETSC_LIBRARIES} 
    DomainCollection
    GMG
    )
if(VTK_FOUND)
    target_link_libraries(test VtkWriter)
endif(VTK_FOUND)
