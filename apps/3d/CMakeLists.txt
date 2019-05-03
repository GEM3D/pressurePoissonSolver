add_executable(steady steady.cpp)
target_link_libraries(steady
    UTILS
    Thunderegg
    tpl
)
