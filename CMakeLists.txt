cmake_minimum_required(VERSION 3.20.3)

project("N_Body_Simulation")
enable_language(C)
enable_language(CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -march=native")


# Get version from Git tags
execute_process(
        COMMAND git describe --tags --always --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
)

# Fallback if not a git repo
if(NOT GIT_VERSION)
    set(GIT_VERSION "Unknown_Version")
endif()

# Get current commit hash
execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
)

# Check for uncommitted changes
execute_process(
        COMMAND git diff-index --quiet HEAD --
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE GIT_DIRTY
        ERROR_QUIET
)

# Set dirty flag
if(GIT_DIRTY)
    set(GIT_DIRTY_FLAG "-dirty")
    message(WARNING "Building with uncommitted changes! Build is NOT reproducible.")

    #    # Save uncommitted changes to a patch file
    #    execute_process(
    #            COMMAND git diff HEAD
    #            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    #            OUTPUT_FILE "${CMAKE_BINARY_DIR}/uncommitted_changes.patch"
    #            ERROR_QUIET
    #    )
    #    message(WARNING "Uncommitted changes saved to: ${CMAKE_BINARY_DIR}/uncommitted_changes.patch")
else()
    set(GIT_DIRTY_FLAG "")
endif()

if(GIT_DIRTY AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "Release builds require committed changes!")
endif()


# Combine to form full version string
set(PROJECT_VERSION "${GIT_VERSION}${GIT_DIRTY_FLAG}")


string(TIMESTAMP BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S")

# Configure build info header
configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/BuildInfo.h.in"
        "${PROJECT_BINARY_DIR}/generated/BuildInfo.h"
        @ONLY
)

#=========================================================================

add_executable(${PROJECT_NAME} src/main.cpp)


set(USE_MPI ON)
set(USE_OPENMP ON)
option(REMOTE_RUN "Configurations for remote run" ON)



# Include dependencies
include(cmake/utils.cmake)
include(cmake/ExternalDependencies.cmake)
include(cmake/PrintDebug.cmake)

add_subdirectory(src/benchmark)


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

