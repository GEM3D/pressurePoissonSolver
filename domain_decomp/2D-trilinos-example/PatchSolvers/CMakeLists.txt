# FFTW
add_library(FftwPatchSolver FftwPatchSolver.cpp)

# Fishpack
enable_language(Fortran)
# Set compiler to convert single precision to double
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
set (CMAKE_Fortran_FLAGS "-r8 -O3 ${CMAKE_Fortran_FLAGS}")
else()
set (CMAKE_Fortran_FLAGS "-fdefault-real-8 -O3 ${CMAKE_Fortran_FLAGS}")
endif()
# Fishpack files
file(GLOB fishpack_SRC "fishpack/*.f")

add_library(FishpackPatchSolver FishpackPatchSolver.cpp ${fishpack_SRC})
target_link_libraries(FftwPatchSolver ${FFTW_LIB})
