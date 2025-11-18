cmake_minimum_required(VERSION 3.20.3)

project("N_Body_Simulation")
set(PROJECT_VERSION "1.0.0")

enable_language(C)
enable_language(CXX)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -march=native")

string(TIMESTAMP CMAKE_TIMESTAMP "%Y-%m-%d %H:%M:%S")

# Configure build info header
configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/BuildInfo.h.in"
        "${PROJECT_BINARY_DIR}/generated/BuildInfo.h"
        @ONLY
)

#=========================================================================

set(USE_MPI ON)
set(USE_OPENMP ON)
option(REMOTE_RUN "Configurations for remote run" ON)



# Include dependencies
include(cmake/utils.cmake)
include(cmake/ExternalDependencies.cmake)
include(cmake/PrintDebug.cmake)

add_subdirectory(src/benchmark)

add_executable(${PROJECT_NAME} src/main.cpp)

if (REMOTE_RUN)
    add_compile_definitions(HPC_RUN)
endif ()

if(USE_MPI)
    find_and_link_mpi_for_target(${PROJECT_NAME})
endif()

if(USE_OPENMP)
    find_and_link_openmp_for_target(${PROJECT_NAME})
endif ()


target_link_libraries(${PROJECT_NAME} PRIVATE print_debug)
target_link_libraries(${PROJECT_NAME} PRIVATE benchmark)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_BINARY_DIR}/generated")

#================================ SET PRINT_DEBUG_LEVEL per module ==========================================

#------------- CUSTOM PRINT DEBUG LEVELS (for console print) -----------------
#[[
 Available log levels - set PRINT_DEBUG_LEVEL to one of:
    PRINT_DEBUG_LEVEL_NONE   → disable all logs
    PRINT_DEBUG_LEVEL_ERROR  → show only errors
    PRINT_DEBUG_LEVEL_WARN   → show warnings and errors
    PRINT_DEBUG_LEVEL_INFO   → show info, warnings, and errors
]]

#------------------------------------------------
# Global debug level (used as default for all modules)
set(PRINT_DEBUG_LEVEL_GLOBAL ${PRINT_DEBUG_LEVEL_INFO} CACHE STRING "Global debug level for all modules")

# Per-module debug levels (use global if not explicitly set)
set(PRINT_DEBUG_LEVEL_MAIN ${PRINT_DEBUG_LEVEL_GLOBAL} CACHE STRING "Debug level for main")
set(PRINT_DEBUG_LEVEL_BENCHMARK ${PRINT_DEBUG_LEVEL_GLOBAL} CACHE STRING "Debug level for benchmark")
#------------------------------------------------

define_print_debug_level_for_target(${PROJECT_NAME} ${PRINT_DEBUG_LEVEL_MAIN})
define_print_debug_level_for_target(benchmark ${PRINT_DEBUG_LEVEL_BENCHMARK})

#================================================================================================

