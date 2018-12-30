project(DomainDecomp)
#
# CMake example that uses FIND_PACKAGE(Trilinos ...) to build your C++
# application with Trilinos.  You should know a little bit about CMake
# before reading this example; in particular, you should know how to
# add C++ source files and header files to your project.
#
# NOTE: This helper function assumes no generator expressions are used
#       for the source files
function(target_sources_local target)
  if(POLICY CMP0076)
    # New behavior is available, so just forward to it by ensuring
    # that we have the policy set to request the new behavior, but
    # don't change the policy setting for the calling scope
    cmake_policy(PUSH)
    cmake_policy(SET CMP0076 NEW)
    target_sources(${target} ${ARGN})
    cmake_policy(POP)
    return()
  endif()

  # Must be using CMake 3.12 or earlier, so simulate the new behavior
  unset(_srcList)
  get_target_property(_targetSourceDir ${target} SOURCE_DIR)

  foreach(src ${ARGN})
    if(NOT src STREQUAL "PRIVATE" AND
       NOT src STREQUAL "PUBLIC" AND
       NOT src STREQUAL "INTERFACE" AND
       NOT IS_ABSOLUTE "${src}")
      # Relative path to source, prepend relative to where target was defined
      file(RELATIVE_PATH src "${_targetSourceDir}" "${CMAKE_CURRENT_LIST_DIR}/${src}")
    endif()
    list(APPEND _srcList ${src})
  endforeach()
  target_sources(${target} ${_srcList})
endfunction()
# You need CMake version >= 2.8 to use this example.
cmake_minimum_required(VERSION 2.8)
enable_language(CXX)
enable_language(C)

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.

set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# file, as above.  This line tells the linker to link the executable 
# with the Trilinos libraries, and with the third-party libraries 
# (TPLs) with which Trilinos was built.
#
#target_link_libraries(steady
#    SchurHelper
#    Init
#    MatrixHelper
#    GMG
#    Interpolators
#    PBMatrix
#    PolyPrecs
#  ${MPI_CXX_LIBRARIES} 
#  ${PETSC_LIBRARIES} 
#  ${FFTW_LIBRARIES} 
#    ${Zoltan_LIBRARIES}
#)
#target_link_libraries(steady2d
#    SchurHelper
#    Init
#    MatrixHelper
#    GMG
#    Interpolators
#    ClawWriter
#    PBMatrix
#    PolyPrecs
#  ${MPI_CXX_LIBRARIES} 
#  ${PETSC_LIBRARIES} 
#  ${FFTW_LIBRARIES} 
#    ${Zoltan_LIBRARIES}
#)
#if(AMGX_FOUND)
#add_library(AmgxWrapper AmgxWrapper.cpp)
#target_link_libraries(steady AmgxWrapper)
#target_link_libraries(AmgxWrapper
#   ${AMGX_LIBRARIES}
#)
#endif()

add_subdirectory(tpl)
add_subdirectory(src)
add_subdirectory(apps)
#SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )
#IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
#  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
#    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
#    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
#  )
#ENDIF()
