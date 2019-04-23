add_library(Thunderegg_GMG Helper.cpp Helper2d.cpp )
target_link_libraries(Thunderegg_GMG PRIVATE Thunderegg)
target_link_libraries(Thunderegg_GMG PRIVATE tpl)
