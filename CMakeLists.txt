cmake_minimum_required(VERSION 3.22)

project(pli CXX)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#==============================================================================
# Configuração do projeto
#==============================================================================
# Leitor de PDF
if (NOT DEFINED PDF_VIEWER)
    set(PDF_VIEWER "open")
endif()

# Caminho para a pasta de instalação do LEMON
if (NOT DEFINED LEMON_DIR)
    set(LEMON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lemon/lemon-1.3.1")
endif()

#==============================================================================
# Configuração da compilação
#==============================================================================
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

#==============================================================================
# Dependências externas
#==============================================================================
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Gurobi - Biblioteca de Otimização
#
# Infere o local de instalação a partir da variável de ambiente GUROBI_HOME ou
# da variável GUROBI_DIR do CMake.
#
# Veja o arquivo cmake/FindGUROBI.cmake para mais detalhes.
find_package(GUROBI REQUIRED)

# LEMON - Biblioteca de otimização em grafos
#
# Infere o local de instalação a partir da variável de ambiente LEMON_HOME ou
# da variável LEMON_DIR do CMake. Por padrão, usa a pasta lemon/lemon-1.3.1
# relativo à raíz do projeto. 
#
# Veja o arquivo cmake/FindLEMON.cmake para mais detalhes.
find_package(LEMON REQUIRED)

#==============================================================================
# Targets
#==============================================================================
configure_file(
    include/network_design/config.hpp.in
    include/network_design/config.hpp
    @ONLY)

include_directories(
    include
    include/mylib
    ${GUROBI_INCLUDE_DIRS}
    ${LEMON_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}/include)

#------------------------------------------------------------------------------
# Biblioteca auxiliar de funções úteis para leitura, visualização e
# manipulação de grafos e afins.
#------------------------------------------------------------------------------
add_library(mylib
    src/mylib/myutils.cpp
    src/mylib/mygraphlib.cpp
    src/mylib/mycolor.cpp
    src/mylib/geompack.cpp
    src/mylib/deprecated.cpp)

target_link_libraries(mylib PUBLIC
    ${LEMON_LIBRARY}
    ${GUROBI_LIBRARIES})

#------------------------------------------------------------------------------
# Biblioteca comum com os algoritmos para solução do problema.
#------------------------------------------------------------------------------
add_library(network_design
    src/network_design/defs.cpp
    src/network_design/heuristics.cpp
    src/network_design/mip/base.cpp
    src/network_design/mip/polynomial.cpp
    src/network_design/mip/branch_and_cut.cpp
    src/network_design/utils.cpp)

target_link_libraries(network_design PUBLIC mylib)

#------------------------------------------------------------------------------
# Executável da parte 1 - Formulação compacta.
#------------------------------------------------------------------------------
add_executable(mo420-network-design-1 src/mo420-network-design-1.cpp)
target_link_libraries(mo420-network-design-1 network_design)

#------------------------------------------------------------------------------
# Executável da parte 2 - Branch-and-Cut.
#------------------------------------------------------------------------------
add_executable(mo420-network-design-2 src/mo420-network-design-2.cpp)
target_link_libraries(mo420-network-design-2 network_design)
