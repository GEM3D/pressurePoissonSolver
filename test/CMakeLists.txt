add_executable(test SchurInfo.cpp Domain.cpp GMG.cpp test.cpp Side.cpp Octant.cpp OctTree.cpp
    DomainCollection.cpp Utils.cpp)
target_link_libraries(test
    Thunderegg
    )
if(VTK_FOUND)
    target_link_libraries(test VtkWriter)
endif(VTK_FOUND)
