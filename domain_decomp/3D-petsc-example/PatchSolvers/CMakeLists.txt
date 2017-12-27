# FFTW
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
add_library(FftwPatchSolver FftwPatchSolver.cpp)
target_link_libraries(FftwPatchSolver ${FFTW_LIB})
