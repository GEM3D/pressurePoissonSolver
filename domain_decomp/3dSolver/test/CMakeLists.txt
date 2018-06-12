add_executable(test SchurDomain.cpp Domain.cpp GMG.cpp test.cpp)
target_link_libraries(test
    ${MPI_CXX_LIBRARIES} 
    ${PETSC_LIBRARIES} 
    DomainCollection
    GMGHelper
    )
