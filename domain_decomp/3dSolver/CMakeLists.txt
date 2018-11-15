#
# CMake example that uses FIND_PACKAGE(Trilinos ...) to build your C++
# application with Trilinos.  You should know a little bit about CMake
# before reading this example; in particular, you should know how to
# add C++ source files and header files to your project.
#

# You need CMake version >= 2.8 to use this example.
cmake_minimum_required(VERSION 2.8)
enable_language(CXX)
enable_language(C)

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.

set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Zoltan REQUIRED)
find_package(PETSc REQUIRED)
# If you haven't already set the C++ compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ need
# to edit this line.
project(DomainDecomp)
find_package(MPI REQUIRED)
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
find_package(VTK)
IF(VTK_FOUND)
    SET(CMAKE_CXX_FLAGS "-DHAVE_VTK ${CMAKE_CXX_FLAGS}")
    MESSAGE("VTK found, enabling vtk output")
ELSE()
    MESSAGE("VTK not found, disabling vtk output")
ENDIF()

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
include_directories(".")
include_directories(${FFTW_INCLUDES})
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_PATH})
if(VTK_FOUND)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
endif()
if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
endif()
if(AMGX_FOUND)
    include_directories(${AMGX_INCLUDE_DIRS})
endif()
include_directories(${PETSC_INCLUDES})
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
add_executable(steady steady.cpp)
add_executable(steady2d steady2d.cpp)
add_library(SchurHelper SchurMatrixHelper.cpp)
add_library(MatrixHelper MatrixHelper.cpp)
add_library(Interpolators TriLinInterp.cpp BilinearInterpolator.cpp)
add_library(Init Init.cpp)
add_library(PBMatrix PBMatrix.cpp)
add_library(PolyPrecs PolyChebPrec.cpp)

#
# EDIT HERE: Change trilinosExample to the name of the executable
# file, as above.  This line tells the linker to link the executable 
# with the Trilinos libraries, and with the third-party libraries 
# (TPLs) with which Trilinos was built.
#
target_link_libraries(steady
    SchurHelper
    Init
    MatrixHelper
    GMG
    Interpolators
    MMWriter
    ClawWriter
    PBMatrix
    PolyPrecs
  ${MPI_CXX_LIBRARIES} 
  ${PETSC_LIBRARIES} 
  ${FFTW_LIBRARIES} 
    ${Zoltan_LIBRARIES}
)
target_link_libraries(steady2d
    SchurHelper
    Init
    MatrixHelper
    GMG
    Interpolators
    MMWriter
    ClawWriter
    PBMatrix
    PolyPrecs
  ${MPI_CXX_LIBRARIES} 
  ${PETSC_LIBRARIES} 
  ${FFTW_LIBRARIES} 
    ${Zoltan_LIBRARIES}
)
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

if(VTK_FOUND)
    target_link_libraries(steady VtkWriter)
    target_link_libraries(steady2d VtkWriter)
endif(VTK_FOUND)
if(AMGX_FOUND)
add_library(AmgxWrapper AmgxWrapper.cpp)
target_link_libraries(steady AmgxWrapper)
target_link_libraries(AmgxWrapper
   ${AMGX_LIBRARIES}
)
endif()

add_subdirectory(PatchSolvers)
add_subdirectory(Writers)
add_subdirectory(test)
add_subdirectory(GMG)
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()
