cmake_minimum_required(VERSION 3.27)

project(pli CXX)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.4)
        message(FATAL_ERROR "GCC version must be at least 9.4!")
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
        message(FATAL_ERROR "Clang version must be at least 10.0!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(GUROBI REQUIRED)
find_package(LEMON REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS} ${LEMON_INCLUDE_DIRS} include include/mylib)

add_library(mylib
    src/mylib/myutils.cpp
    src/mylib/mygraphlib.cpp
    src/mylib/mycolor.cpp
    src/mylib/geompack.cpp
    src/mylib/deprecated.cpp)

add_library(network_design
    src/network_design/defs.cpp
    src/network_design/heuristics.cpp
    src/network_design/mip/base.cpp
    src/network_design/mip/polynomial.cpp
    src/network_design/mip/branch_and_cut.cpp
    src/network_design/utils.cpp)

target_link_libraries(network_design PUBLIC
    mylib
    ${LEMON_LIBRARY}
    ${GUROBI_CXX_LIBRARY} ${GUROBI_C_LIBRARY}
)

add_executable(mo420-network-design-1 src/mo420-network-design-1.cpp)
target_link_libraries(mo420-network-design-1 network_design)

add_executable(mo420-network-design-2 src/mo420-network-design-2.cpp)
target_link_libraries(mo420-network-design-2 network_design)
