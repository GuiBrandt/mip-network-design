find_path(LEMON_INCLUDE_DIRS
    NAMES lemon/config.h
    HINTS ${LEMON_DIR} $ENV{LEMON_HOME}
    PATH_SUFFIXES include
    DOC "LEMON include path")
    
find_library(LEMON_LIBRARY
    NAMES emon
    HINTS ${LEMON_DIR} $ENV{LEMON_HOME}
    PATH_SUFFIXES lib
    DOC "LEMON library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LEMON DEFAULT_MSG LEMON_LIBRARY)
