project(Thunderegg_Lib)

find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)
find_package(FFTW REQUIRED)
find_package(Zoltan REQUIRED)
find_package(CUDA)
find_package(p4est)
if(p4est_FOUND)
    find_package(sc REQUIRED)
endif()

#determine sources first
list(APPEND Thunderegg_SRCS Thunderegg/SchurMatrixHelper.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/SchurMatrixHelper2d.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/MatrixHelper.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/MatrixHelper2d.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/TriLinInterp.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/BilinearInterpolator.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/PBMatrix.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/PolyChebPrec.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/SchurHelper.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/DomainCollection.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/BalancedLevelsGenerator.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/ThundereggDCG.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/Domain.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/SchurDomain.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/Vector.cpp)

list(APPEND Thunderegg_SRCS Thunderegg/GMG/Helper.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/GMG/Helper2d.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/GMG/CycleFactory2d.cpp)

list(APPEND Thunderegg_SRCS Thunderegg/PatchSolvers/DftPatchSolver.cpp)
list(APPEND Thunderegg_SRCS Thunderegg/PatchSolvers/FftwPatchSolver.cpp)

if(p4est_FOUND)
    list(APPEND Thunderegg_SRCS Thunderegg/p4estBLG.cpp)
    list(APPEND Thunderegg_SRCS Thunderegg/P4estDCG.cpp)
endif()

add_library(Thunderegg SHARED ${Thunderegg_SRCS})

# set include directories
target_include_directories(Thunderegg PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
list(APPEND Thunderegg_Includes ${MPI_CXX_INCLUDE_PATH})
list(APPEND Thunderegg_Includes ${PETSC_INCLUDES})
list(APPEND Thunderegg_Includes ${FFTW_INCLUDES})
list(APPEND Thunderegg_Includes ${Zoltan_INCLUDES})
if(p4est_FOUND)
    list(APPEND Thunderegg_Includes ${p4est_INCLUDES})
    list(APPEND Thunderegg_Includes ${sc_INCLUDES})
endif()
target_include_directories(Thunderegg PUBLIC ${Thunderegg_Includes})

# set link 
list(APPEND Thunderegg_Libs ${MPI_CXX_LIBRARIES})
list(APPEND Thunderegg_Libs ${PETSC_LIBRARIES}) 
list(APPEND Thunderegg_Libs ${FFTW_LIBRARIES}) 
list(APPEND Thunderegg_Libs ${Zoltan_LIBRARIES})
if(p4est_FOUND)
    list(APPEND Thunderegg_Libs ${p4est_LIBRARIES})
    list(APPEND Thunderegg_Libs ${sc_LIBRARIES})
endif()
target_link_libraries(Thunderegg PUBLIC ${Thunderegg_Libs})


install(TARGETS Thunderegg EXPORT ThundereggConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})


install(DIRECTORY Thunderegg/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Thunderegg
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

install(EXPORT ThundereggConfig DESTINATION share/Thunderegg/cmake)

foreach(LIB ${Thunderegg_Libs})
    set(Thunderegg_Lib_Flags "${Thunderegg_Lib_Flags} -l${LIB}")
endforeach()

foreach(LIB ${Thunderegg_Includes})
    set(Thunderegg_Include_Flags "${Thunderegg_Include_Flags} -I${LIB}")
endforeach()

CONFIGURE_FILE("Thunderegg.pc.in" "Thunderegg.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Thunderegg.pc" 
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

