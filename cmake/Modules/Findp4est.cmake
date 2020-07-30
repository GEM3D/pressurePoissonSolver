# - Try to find p4est
#

find_path (p4est_DIR include/p4est.h HINTS p4est_DIR ENV p4est_DIR CPATH)

  SET(p4est_INCLUDES ${p4est_DIR})
  find_path (p4est_INCLUDE_DIR p4est.h HINTS "${p4est_DIR}/include" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND p4est_INCLUDES ${p4est_INCLUDE_DIR})
  find_library(p4est_LIBRARIES p4est PATHS "${p4est_DIR}/lib" ${p4est_DIR})

  #get linker flags from config header
  file(READ ${p4est_INCLUDE_DIR}/p4est_config.h p4est_config_h)
  string(REGEX MATCH "#define P4EST_LIBS [^\n]*" p4est_extra_lib_string ${p4est_config_h})
  string(REGEX MATCH "\".*\"" p4est_extra_lib_string ${p4est_extra_lib_string})
  message("-- found p4est linker flags: ${p4est_extra_lib_string}")

  #get library dir from linker flags
  string(REGEX MATCHALL "-L[^ ]*" p4est_tmp ${p4est_extra_lib_string})
  foreach(x ${p4est_tmp})
    string(SUBSTRING ${x} 2 -1 y)
    list(APPEND p4est_extra_lib_dirs ${y})
  endforeach(x)
  message("-- found p4est dependant library_dirs: ${p4est_extra_lib_dirs}")

  #get libraries from linker flags
  string(REGEX MATCHALL "-l[^ ]*" p4est_tmp ${p4est_extra_lib_string})
  foreach(x ${p4est_tmp})
    string(SUBSTRING ${x} 2 -1 y)
    list(APPEND p4est_extra_libs ${y})
  endforeach(x)
  message("-- found p4est dependant libraries: ${p4est_extra_libs}")

  foreach(lib ${p4est_extra_libs})
    find_library(extra_lib ${lib} ${p4est_extra_lib_dirs})
    list(APPEND p4est_LIBRARIES ${extra_lib})
    unset(extra_lib CACHE)
  endforeach(lib)
  
  message("-- p4est libraries: ${p4est_LIBRARIES}")

  list(APPEND p4est_LIBRARIES ${p4est_extra_lib})

IF(EXISTS ${p4est_DIR}/include/p4est.h)
  SET(p4est_FOUND YES)
ELSE(EXISTS ${p4est_DIR}/include/p4est.h)
  SET(p4est_FOUND NO)
ENDIF(EXISTS ${p4est_DIR}/include/p4est.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p4est DEFAULT_MSG p4est_LIBRARIES p4est_INCLUDES)
