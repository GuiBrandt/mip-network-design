find_path(GUROBI_INCLUDE_DIRS
    NAMES gurobi_c++.h
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES include
    DOC "GUROBI include path")

find_library(GUROBI_CXX_LIBRARY
    NAMES gurobi_c++
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib
    DOC "GUROBI C++ library")

find_library(GUROBI_C_LIBRARY
    NAMES gurobi100
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib
    DOC "GUROBI C library")

set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_CXX_LIBRARY)
