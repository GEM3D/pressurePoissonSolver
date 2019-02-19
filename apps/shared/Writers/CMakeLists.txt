target_sources_local(UTILS PRIVATE ClawWriter.cpp ClawWriter.h)
if(VTK_FOUND)
    target_sources_local(UTILS PRIVATE VtkWriter.cpp VtkWriter2d.cpp VtkWriter.h
        VtkWriter2d.h)
    link_directories(${VTK_LIBRARY_DIRS})
    target_include_directories(UTILS PUBLIC ${VTK_INCLUDE_DIRS})
    target_link_libraries(UTILS ${VTK_LIBRARIES})
endif()
