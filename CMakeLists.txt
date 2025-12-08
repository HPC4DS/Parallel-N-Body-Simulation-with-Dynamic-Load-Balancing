cmake_minimum_required(VERSION 3.20.3)

project("N_Body_Simulation")
set(PROJECT_VERSION "1.0.0")

enable_language(C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
#set(CMAKE_C_EXTENSIONS OFF)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)


string(TIMESTAMP CMAKE_TIMESTAMP "%Y-%m-%d %H:%M:%S")

# Configure build info header
configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/BuildInfo.h.in"
        "${PROJECT_BINARY_DIR}/generated/BuildInfo.h"
        @ONLY
)

#======================== OPTIONS ==============================
set(USE_MPI ON)
set(USE_OPENMP ON)
option(REMOTE_RUN "Configurations for remote run" ON)
option(BUILD_TESTING "Enable testing" OFF)

add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

if (REMOTE_RUN)
    add_compile_definitions(HPC_RUN)
endif ()

#======================== SET PRINT_DEBUG_LEVEL ==========================

#------------- CUSTOM PRINT DEBUG LEVELS  -----------------
#[[
 Available log levels - set PRINT_DEBUG_LEVEL to one of:
    PRINT_DEBUG_LEVEL_NONE   → disable all logs
    PRINT_DEBUG_LEVEL_ERROR  → show only errors
    PRINT_DEBUG_LEVEL_WARN   → show warnings and errors
    PRINT_DEBUG_LEVEL_INFO   → show info, warnings, and errors
]]
#----------------------------------------------------------
# Global debug level (used as default for all modules)
set(PRINT_DEBUG_LEVEL_GLOBAL ${PRINT_DEBUG_LEVEL_INFO} CACHE STRING "Global debug level for all modules")

#----------------------------------------------------------
# Available print debug modules:
# APP_PRINT_DEBUG_LEVEL
# BENCHMARK_PRINT_DEBUG_LEVEL
#----------------------------------------------------------
#=========================================================================


#======================== Include dependencies ===========================

include(cmake/cmake_utils.cmake)
include(cmake/ExternalDependencies.cmake)
include(cmake/PrintDebug.cmake)
include(cmake/utils.cmake)


add_subdirectory(src)
add_subdirectory(benchmarks/src)
add_subdirectory(benchmarks/microbenchmarks)
add_subdirectory(app)
add_subdirectory(sandbox)

if(BUILD_TESTING)
    include(cmake/fetch_googletest.cmake)
    add_subdirectory(tests)
endif()

#==============================================
# UTILITIES TO COPY IN OTHER CMAKE FILES
#==============================================

#if(USE_MPI)
#    find_and_link_mpi_for_target(${PROJECT_NAME})
#endif()
#
#if(USE_OPENMP)
#    find_and_link_openmp_for_target(${PROJECT_NAME})
#endif ()


# target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_BINARY_DIR}/generated")


# ------
# This shouldn't be needed anywhere, but just in case...
# ------
#include(${CMAKE_SOURCE_DIR}/cmake/PrintDebug.cmake)
#include(${CMAKE_SOURCE_DIR}/cmake/ExternalDependencies.cmake)
#include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
