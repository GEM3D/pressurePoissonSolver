project(Thunderegg_Lib)

find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)
find_package(FFTW REQUIRED)
find_package(Zoltan REQUIRED)
find_package(CUDA)
find_package(p4est)

#determine sources first
list(APPEND Thunderegg_SRCS SchurMatrixHelper.cpp SchurMatrixHelper2d.cpp)
list(APPEND Thunderegg_SRCS MatrixHelper.cpp MatrixHelper2d.cpp)
list(APPEND Thunderegg_SRCS TriLinInterp.cpp BilinearInterpolator.cpp)
list(APPEND Thunderegg_SRCS PBMatrix.cpp)
list(APPEND Thunderegg_SRCS PolyChebPrec.cpp)
if(p4est_FOUND)
    list(APPEND Thunderegg_SRCS p4estBLG.cpp)
endif()

add_library(Thunderegg ${Thunderegg_SRCS})

# set include directories
target_include_directories(Thunderegg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Thunderegg PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_include_directories(Thunderegg PUBLIC ${PETSC_INCLUDES})
target_include_directories(Thunderegg PUBLIC ${FFTW_INCLUDES})
target_include_directories(Thunderegg PUBLIC ${Zoltan_INCLUDES})
if(p4est_FOUND)
    target_include_directories(Thunderegg PUBLIC ${p4est_INCLUDES})
endif()

# set link 
target_link_libraries(Thunderegg PUBLIC tpl)
target_link_libraries(Thunderegg PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(Thunderegg PUBLIC ${PETSC_LIBRARIES}) 
target_link_libraries(Thunderegg PUBLIC ${FFTW_LIBRARIES}) 
target_link_libraries(Thunderegg PUBLIC ${Zoltan_LIBRARIES})
if(p4est_FOUND)
    target_link_libraries(Thunderegg PUBLIC ${p4est_LIBRARIES})
endif()
target_link_libraries(Thunderegg PUBLIC Thunderegg_GMG)

add_subdirectory(PatchSolvers)
add_subdirectory(GMG)
