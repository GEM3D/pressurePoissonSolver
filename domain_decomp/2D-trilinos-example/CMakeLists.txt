#
# CMake example that uses FIND_PACKAGE(Trilinos ...) to build your C++
# application with Trilinos.  You should know a little bit about CMake
# before reading this example; in particular, you should know how to
# add C++ source files and header files to your project.
#
enable_language(Fortran)

# You need CMake version >= 2.8 to use this example.
cmake_minimum_required(VERSION 2.8)

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.
find_package(Trilinos PATHS $ENV{TRILINOS_DIR}/lib/cmake/Trilinos $ENV{TRILINOS_DIR})

# If FIND_PACKAGE successfully found your Trilinos install, it will
# set the Boolean flag Trilinos_FOUND.  The following IF statement
# fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
# found, it prints out the values of some Trilinos configuration
# details.  You may find them useful for building your application
# that uses Trilinos.
IF(Trilinos_FOUND)
  MESSAGE("Found Trilinos!")
ELSE()
    MESSAGE(FATAL_ERROR "Could not find Trilinos!\n Is TRILINOS_DIR set?")
ENDIF()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# If you haven't already set the C++ compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ need
# to edit this line.
project(DomainDecomp)
find_package(FFTW REQUIRED)
find_package(MPI REQUIRED)
set (CMAKE_CXX_FLAGS "--std=c++11 -Wall --pedantic -O3 ${CMAKE_CXX_FLAGS}")
IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
ENDIF()

# If you haven't already set the C compiler, use the same compiler
# that was used to build Trilinos on your machine.  You do _not_ 
# need to edit this line.
IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILE ${Trilinos_C_COMPILER})
ENDIF()

set (CMAKE_Fortran_FLAGS "-fdefault-real-8 ${CMAKE_Fortran_FLAGS}")

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  MESSAGE("Using ccache.")
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)
#
# EDIT HERE: Replace TrilinosExample with your project's name.
#

# 
# EDIT HERE: 
#
# Source file(s) in your project.  Here, we include one source file, 
# main.cpp, which contains the main() function.  You may also want to 
# set HEADERS to the list of header files in your project, if there 
# are any.
#
set(SOURCES 
  heat.cpp
  )

# 
# EDIT HERE:
#
# Add to the header include path the current directory, Trilinos'
# include directory, and the include directories of the third-party
# libraries (TPLs) with which Trilinos was built.
#
include_directories(".")
include_directories(${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${FFTW_INCLUDES})
include_directories(${MPI_CXX_INCLUDE_PATH})

# Add to the library path Trilinos' library path, and the library
# paths of the third-party libraries (TPLs) with which Trilinos was
# built.  
link_directories(${Trilinos_LIBRARY_DIRS})
link_directories(${Trilinos_TPL_LIBRARY_DIRS})

# 
# EDIT HERE: Change trilinosExample to the name of the executable
# file that you want to build.  This line tells CMake which source
# file(s) to use to build the executable.
#
add_executable(heat heat.cpp)
add_library(DomainCollection DomainCollection.cpp)
add_library(DomainSignatureCollection DomainSignatureCollection.cpp)
add_library(Domain Domain.cpp)
add_library(RBMatrix RBMatrix.cpp)
add_library(Factory Factory.cpp)
add_library(FftwSolver FftwSolver.cpp)
file(GLOB fishpack_SRC "fishpack/*.f")
#add_library(FishpackSolver FishpackSolver.cpp)
#add_library(Fishpack ${fishpack_SRC})

#
# EDIT HERE: Change trilinosExample to the name of the executable
# file, as above.  This line tells the linker to link the executable 
# with the Trilinos libraries, and with the third-party libraries 
# (TPLs) with which Trilinos was built.
#
target_link_libraries(heat
    DomainCollection
    DomainSignatureCollection
    RBMatrix
    Domain
    Factory
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
  ${FFTW_LIB}
  ${MPI_CXX_LIBRARIES} 
)
target_link_libraries(Factory
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
  ${MPI_CXX_LIBRARIES} 
)
target_link_libraries(RBMatrix
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
  ${FFTW_LIB}
  ${MPI_CXX_LIBRARIES} 
)
target_link_libraries(DomainCollection
    RBMatrix 
    Domain
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
  ${FFTW_LIB}
  ${MPI_CXX_LIBRARIES} 
)
target_link_libraries(Domain
    FftwSolver
    #    FishpackSolver
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
  ${FFTW_LIB}
  ${MPI_CXX_LIBRARIES} 
)
#target_link_libraries(FishpackSolver
#    Domain
#    Fishpack
#)
target_link_libraries(FftwSolver
    Domain
)
target_link_libraries(DomainSignatureCollection
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES} 
)
