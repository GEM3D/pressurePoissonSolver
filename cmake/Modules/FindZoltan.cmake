# - Try to find Zoltan
#

find_path (ZOLTAN_DIR include/zoltan.h HINTS ZOLTAN_DIR ENV ZOLTAN_DIR CPATH)

  SET(Zoltan_INCLUDES ${ZOLTAN_DIR})
  find_path (Zoltan_INCLUDE_DIR zoltan.h HINTS "${ZOLTAN_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND Zoltan_INCLUDES ${Zoltan_INCLUDE_DIR})
  find_library(Zoltan_LIBRARIES zoltan PATHS "${ZOLTAN_DIR}/lib" ${ZOLTAN_DIR})
IF(EXISTS ${ZOLTAN_DIR}/include/zoltan.h)
  SET(Zoltan_FOUND YES)
ELSE(EXISTS ${ZOLTAN_DIR}/include/zoltan.h)
  SET(Zoltan_FOUND NO)
ENDIF(EXISTS ${ZOLTAN_DIR}/include/zoltan.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan DEFAULT_MSG Zoltan_LIBRARIES Zoltan_INCLUDES)
