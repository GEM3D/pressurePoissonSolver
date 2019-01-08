project(DomainDecomp)
add_executable(steady2d steady.cpp)
target_link_libraries(steady2d
    PPS
    UTILS
)
