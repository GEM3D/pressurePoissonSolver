# - Try to find p4est
#

find_path (p4est_DIR include/p4est.h HINTS p4est_DIR ENV p4est_DIR CPATH)

  SET(p4est_INCLUDES ${p4est_DIR})
  find_path (p4est_INCLUDE_DIR p4est.h HINTS "${p4est_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND p4est_INCLUDES ${p4est_INCLUDE_DIR})
  find_library(p4est_LIBRARIES p4est PATHS "${p4est_DIR}/lib" ${p4est_DIR})
IF(EXISTS ${p4est_DIR}/include/zoltan.h)
  SET(p4est_FOUND YES)
ELSE(EXISTS ${p4est_DIR}/include/zoltan.h)
  SET(p4est_FOUND NO)
ENDIF(EXISTS ${p4est_DIR}/include/zoltan.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p4est DEFAULT_MSG p4est_LIBRARIES p4est_INCLUDES)
