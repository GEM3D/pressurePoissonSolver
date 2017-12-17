# - Try to find Zoltan
#

find_path (Zoltan_DIR include/zoltan_cpp.h HINTS ENV Zoltan_DIR)

IF(EXISTS ${Zoltan_DIR}/include/zoltan_cpp.h)
  SET(Zoltan_FOUND YES)
  SET(Zoltan_INCLUDES ${Zoltan_DIR})
  find_path (Zoltan_INCLUDE_DIR zoltan_cpp.h HINTS "${Zoltan_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND Zoltan_INCLUDES ${Zoltan_INCLUDE_DIR})
  find_library(Zoltan_LIBRARIES zoltan PATHS "${Zoltan_DIR}/lib" ${Zoltan_DIR})
ELSE(EXISTS ${Zoltan_DIR}/include/zoltan_cpp.h)
  SET(Zoltan_FOUND NO)
ENDIF(EXISTS ${Zoltan_DIR}/include/zoltan_cpp.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan DEFAULT_MSG Zoltan_LIBRARIES Zoltan_INCLUDES)
