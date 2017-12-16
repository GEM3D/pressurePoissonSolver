# FFTW
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
add_library(FftwPatchSolver FftwPatchSolver.cpp)

# Fishpack
enable_language(Fortran)
# Set compiler to convert single precision to double
set (CMAKE_Fortran_FLAGS "-fdefault-real-8 -O3 ${CMAKE_Fortran_FLAGS}")
# Fishpack files
file(GLOB fishpack_SRC "fishpack/*.f")

add_library(FishpackPatchSolver FishpackPatchSolver.cpp ${fishpack_SRC})
target_link_libraries(FftwPatchSolver ${FFTW_LIB})
