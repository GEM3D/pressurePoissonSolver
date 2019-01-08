project(DomainDecomp)

find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)
set (CMAKE_CXX_FLAGS "-Wall --std=c++11 -O3 -march=native ${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG  "-Wall --std=c++11 -g -O3  -DDD_DEBUG")
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# If you haven't already set the C compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ 
# need to edit this line.
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILE ${MPI_C_COMPILER})

find_package(VTK)
IF(VTK_FOUND)
    SET(CMAKE_CXX_FLAGS "-DHAVE_VTK ${CMAKE_CXX_FLAGS}")
    MESSAGE("VTK found, enabling vtk output")
ELSE()
    MESSAGE("VTK not found, disabling vtk output")
ENDIF()

include_directories(".")
include_directories(${MPI_CXX_INCLUDE_PATH})
if(VTK_FOUND)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
endif()
include_directories(${PETSC_INCLUDES})

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
add_library(UTILS "")
target_sources_local(UTILS PRIVATE Init.cpp)
target_link_libraries(UTILS PPS)
target_include_directories(UTILS PUBLIC ".")

include(Writers/CMakeLists.txt)
