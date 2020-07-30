# - Try to find Zoltan
#

find_path (Zoltan_DIR include/zoltan.h HINTS ENV Zoltan_DIR)

#find libs
find_library(
  Zoltan_LIB
  NAMES "zoltan"
  PATHS ${Zoltan_DIR}
  PATH_SUFFIXES "lib" "lib64"
)

#find includes
find_path(
  Zoltan_INCLUDES
  NAMES "zoltan.h"
  PATHS ${Zoltan_DIR}
  PATH_SUFFIXES "include"
)


set(Zoltan_LIBRARIES ${Zoltan_LIB} ${ZoltanF_LIB})

if(ZoltanL_LIB)
  set(Zoltan_LIBRARIES ${Zoltan_LIBRARIES} ${ZoltanL_LIB})
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan DEFAULT_MSG
                                  Zoltan_INCLUDES Zoltan_LIBRARIES)

mark_as_advanced(Zoltan_INCLUDES Zoltan_LIBRARIES Zoltan_LIB )

