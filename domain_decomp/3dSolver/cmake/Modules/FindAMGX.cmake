# - Try to find AMGX
#

find_path (AMGX_DIR include/amgx_c.h HINTS ENV AMGX_DIR)

IF(EXISTS ${AMGX_DIR}/include/amgx_c.h)
  SET(AMGX_FOUND YES)
  SET(AMGX_INCLUDES ${AMGX_DIR})
  find_path (AMGX_INCLUDE_DIR amgx_c.h HINTS "${AMGX_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND AMGX_INCLUDES ${AMGX_INCLUDE_DIR})
  FILE(GLOB AMGX_LIBRARIES RELATIVE "${AMGX_DIR}/lib" "${AMGX_DIR}/lib/libamgxsh.so")
ELSE(EXISTS ${AMGX_DIR}/include/amgx_c.h)
  SET(AMGX_FOUND NO)
ENDIF(EXISTS ${AMGX_DIR}/include/amgx_c.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMGX DEFAULT_MSG AMGX_LIBRARIES AMGX_INCLUDES)
