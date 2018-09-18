# FFTW
add_library(FftwPatchSolver FftwPatchSolver.cpp DftPatchSolver.cpp)
target_link_libraries(FftwPatchSolver ${FFTW_LIB})
