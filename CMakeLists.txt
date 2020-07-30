cmake_minimum_required(VERSION 2.8)
project(Thunderegg)
enable_language(CXX)
enable_language(C)

#set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (CMAKE_CXX_FLAGS "-Wall --std=c++11 -O3 -march=native ${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG  "-Wall --std=c++11 -g -O3 -march=native -DDD_DEBUG")

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  MESSAGE("-- Using ccache.")
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

add_subdirectory(tpl)
add_subdirectory(src)
add_subdirectory(apps)
if(CMAKE_BUILD_TYPE EQUAL Debug)
add_subdirectory(test)
endif(CMAKE_BUILD_TYPE EQUAL Debug)
