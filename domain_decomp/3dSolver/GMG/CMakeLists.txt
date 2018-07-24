add_library(GMG Helper.cpp AvgRstr.cpp DrctIntp.cpp TriLinIntp.cpp Cycle.cpp
    InterLevelComm.cpp)
target_link_libraries(GMG
    DomainCollection
)

