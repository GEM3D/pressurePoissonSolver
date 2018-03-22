# FFTW
add_library(FftwPatchSolver FftwPatchSolver.cpp)
target_link_libraries(FftwPatchSolver ${FFTW_LIB})
