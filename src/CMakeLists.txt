project(DomainDecomp)

find_package(Zoltan REQUIRED)
find_package(PETSc REQUIRED)
find_package(p4est REQUIRED)
find_package(MPI REQUIRED)
# If you haven't already set the C++ compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ need
# to edit this line.
set (CMAKE_CXX_FLAGS "-Wall --std=c++11 -O3 -march=native ${CMAKE_CXX_FLAGS}")
if(NUM_GPU)
    set (CMAKE_CXX_FLAGS "-DNUM_GPU=${NUM_GPU} ${CMAKE_CXX_FLAGS}")
endif()
set (CMAKE_CXX_FLAGS_DEBUG  "-Wall --std=c++11 -g -O3  -DDD_DEBUG")
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

find_package(CUDA)
find_package(AMGX)
if(AMGX_FOUND)
    SET(CMAKE_CXX_FLAGS "-DENABLE_AMGX ${CMAKE_CXX_FLAGS}")
endif()
# If you haven't already set the C compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ 
# need to edit this line.
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILE ${MPI_C_COMPILER})


find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  MESSAGE("Using ccache.")
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

# 
# EDIT HERE:
#
# Add to the header include path the current directory, Trilinos'
# include directory, and the include directories of the third-party
# libraries (TPLs) with which Trilinos was built.
#
find_package(FFTW REQUIRED)
add_library(Thunderegg "")
include_directories(".")
include_directories(${FFTW_INCLUDES})
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_PATH})
if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
endif()
if(AMGX_FOUND)
    include_directories(${AMGX_INCLUDE_DIRS})
endif()
target_include_directories(Thunderegg PUBLIC ${PETSC_INCLUDES}
    ${FFTW_INCLUDES} 
    ${Zoltan_INCLUDES}
    )
include_directories(${Zoltan_INCLUDES})

# Add to the library path Trilinos' library path, and the library
# paths of the third-party libraries (TPLs) with which Trilinos was
# built.  
if(AMGX_FOUND)
link_directories("${AMGX_DIR}/lib")
endif()
link_directories("${Zoltan_DIR}/lib")

# 
# EDIT HERE: Change trilinosExample to the name of the executable
# file that you want to build.  This line tells CMake which source
# file(s) to use to build the executable.
#
target_include_directories(Thunderegg PUBLIC ".")
target_sources_local(Thunderegg PRIVATE SchurMatrixHelper.cpp SchurMatrixHelper2d.cpp)
target_sources_local(Thunderegg PRIVATE MatrixHelper.cpp MatrixHelper2d.cpp)
target_sources_local(Thunderegg PRIVATE TriLinInterp.cpp BilinearInterpolator.cpp)
target_sources_local(Thunderegg PRIVATE PBMatrix.cpp)
target_sources_local(Thunderegg PRIVATE PolyChebPrec.cpp)
target_link_libraries(
    Thunderegg PUBLIC
    tpl
    ${MPI_CXX_LIBRARIES} 
    ${PETSC_LIBRARIES} 
    ${FFTW_LIBRARIES} 
    ${Zoltan_LIBRARIES}
)
#
# EDIT HERE: Change trilinosExample to the name of the executable
# file, as above.  This line tells the linker to link the executable 
# with the Trilinos libraries, and with the third-party libraries 
# (TPLs) with which Trilinos was built.
#
if(ENABLE_MUELU)
    find_library(MUELU_LIB MueLuWrapper)
    include_directories("./lib/")
    target_link_libraries(steady
        ${MUELU_LIB}
    )
    SET(CMAKE_CXX_FLAGS "-DENABLE_MUELU ${CMAKE_CXX_FLAGS}")
endif()
if(ENABLE_MUELU_CUDA)
    target_link_libraries(steady
        ${CUDA_LIBRARIES}
        ${CUDA_CUBLAS_LIBRARIES}
    )
    SET(CMAKE_CXX_FLAGS "-DENABLE_MUELU_CUDA ${CMAKE_CXX_FLAGS}")
endif()

if(AMGX_FOUND)
add_library(AmgxWrapper AmgxWrapper.cpp)
target_link_libraries(steady AmgxWrapper)
target_link_libraries(AmgxWrapper
   ${AMGX_LIBRARIES}
)
endif()


add_subdirectory(PatchSolvers)
include(GMG/CMakeLists.txt)
