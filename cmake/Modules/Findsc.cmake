# - Try to find sc
#

find_path (sc_DIR include/sc.h HINTS p4est_DIR sc_DIR ENV p4est_DIR ENV sc_DIR CPATH)

  SET(sc_INCLUDES ${sc_DIR})
  find_path (sc_INCLUDE_DIR sc.h HINTS "${sc_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND sc_INCLUDES ${sc_INCLUDE_DIR})
  find_library(sc_LIBRARIES sc PATHS "${sc_DIR}/lib" ${sc_DIR})
IF(EXISTS ${sc_DIR}/include/sc.h)
  SET(sc_FOUND YES)
ELSE(EXISTS ${sc_DIR}/include/sc.h)
  SET(sc_FOUND NO)
ENDIF(EXISTS ${sc_DIR}/include/sc.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sc DEFAULT_MSG sc_LIBRARIES sc_INCLUDES)
