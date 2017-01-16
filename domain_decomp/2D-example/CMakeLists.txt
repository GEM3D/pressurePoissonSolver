cmake_minimum_required(VERSION 2.8)

# Project's Name and Lanuages
project(2D_Example CXX C)
set (CMAKE_CXX_FLAGS "--std=c++11 -O2 ${CMAKE_CXX_FLAGS}")

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(Boost REQUIRED)
find_package(FFTW REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDES})
# includes and libraries
#link_directories(${FFTW_LIBRARIES})
add_executable(heat heat.cpp)
add_library(Domain Domain.cpp)
TARGET_LINK_LIBRARIES(heat Domain)
TARGET_LINK_LIBRARIES(Domain ${FFTW_LIB})
