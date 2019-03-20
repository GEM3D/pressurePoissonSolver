project(DomainDecomp)

find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)
set (CMAKE_CXX_FLAGS "-Wall --std=c++11 -O3 -march=native ${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS "-Wall --std=c++11 -O3 -march=native ${CMAKE_CXX_FLAGS}" PARENT_SCOPE)
set (CMAKE_CXX_FLAGS_DEBUG  "-Wall --std=c++11 -g -O3  -DDD_DEBUG")
set (CMAKE_CXX_FLAGS_DEBUG  "-Wall --std=c++11 -g -O3  -DDD_DEBUG" PARENT_SCOPE)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# If you haven't already set the C compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ 
# need to edit this line.
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILE ${MPI_C_COMPILER})

find_package(VTK)
IF(VTK_FOUND)
    SET(CMAKE_CXX_FLAGS "-DHAVE_VTK ${CMAKE_CXX_FLAGS}" PARENT_SCOPE)
    MESSAGE("VTK found, enabling vtk output")
ELSE()
    MESSAGE("VTK not found, disabling vtk output")
ENDIF()


add_library(UTILS "")
target_sources_local(UTILS PRIVATE Init.cpp Init.h)
target_link_libraries(UTILS Thunderegg)
target_include_directories(UTILS PUBLIC ".")

find_package(p4est REQUIRED)
target_sources_local(UTILS PRIVATE TreeToP4est.cpp TreeToP4est.h)
target_include_directories(UTILS PUBLIC ${p4est_INCLUDES})
target_link_libraries(UTILS ${p4est_LIBRARIES})

include(Writers/CMakeLists.txt)
