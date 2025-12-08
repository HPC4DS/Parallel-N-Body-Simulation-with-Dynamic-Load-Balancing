//
// Created by Matteo Ranzi on 19/11/25.
//

#ifndef N_BODY_SIMULATION_BUILD_UTILS_H
#define N_BODY_SIMULATION_BUILD_UTILS_H

#include "BuildInfo.h"
#include "debug/unique_print_debug.h"

#ifdef __cplusplus
extern "C" {
#endif


static inline void print_build_info(const int rank) {
    UNIQUE_PRINT(rank, "=============== Build Information ================\n");
#ifdef __cplusplus
    UNIQUE_PRINT(rank, "CXX Standard: %ld\n", __cplusplus);
#else
    UNIQUE_PRINT(rank, "C Standard: %ld\n", __STDC_VERSION__);
#endif
    UNIQUE_PRINT(rank, "Build Version: %s\n", BUILD_VERSION);
    UNIQUE_PRINT(rank, "Build Type: %s\n", BUILD_TYPE);
    UNIQUE_PRINT(rank, "Compiler: %s\n", CXX_COMPILER);
    UNIQUE_PRINT(rank, "CXX Flags: %s\n", CXX_FLAGS);
    UNIQUE_PRINT(rank, "CXX Release Flags: %s\n", CXX_FLAGS_RELEASE);
    UNIQUE_PRINT(rank, "CXX Debug Flags: %s\n", CXX_FLAGS_DEBUG);
    UNIQUE_PRINT(rank, "C Flags: %s\n", C_FLAGS);
    UNIQUE_PRINT(rank, "C Release Flags: %s\n", C_FLAGS_RELEASE);
    UNIQUE_PRINT(rank, "C Debug Flags: %s\n", C_FLAGS_DEBUG);
    UNIQUE_PRINT(rank, "Build Host: %s\n", BUILD_HOST);
    UNIQUE_PRINT(rank, "CMake Timestamp: %s\n", CMAKE_TIMESTAMP);
    UNIQUE_PRINT(rank, "System: %s (%s)\n", SYSTEM_NAME, SYSTEM_PROCESSOR);
    UNIQUE_PRINT(rank, "==================================================\n\n");

}

static inline void build_info_to_string(char* buffer, const size_t buffer_size) {
    int offset = 0;
#ifdef __cplusplus
    offset += snprintf(buffer + offset, buffer_size - offset, "CXX Standard: %ld\n", __cplusplus);
#else
    offset += snprintf(buffer + offset, buffer_size - offset, "C Standard: %ld\n", __STDC_VERSION__);
#endif
    offset += snprintf(buffer + offset, buffer_size - offset, "Build Version: %s\n", BUILD_VERSION);
    offset += snprintf(buffer + offset, buffer_size - offset, "Build Type: %s\n", BUILD_TYPE);
    offset += snprintf(buffer + offset, buffer_size - offset, "Compiler: %s\n", CXX_COMPILER);
    offset += snprintf(buffer + offset, buffer_size - offset, "CXX Flags: %s\n", CXX_FLAGS);
    offset += snprintf(buffer + offset, buffer_size - offset, "CXX Release Flags: %s\n", CXX_FLAGS_RELEASE);
    offset += snprintf(buffer + offset, buffer_size - offset, "CXX Debug Flags: %s\n", CXX_FLAGS_DEBUG);
    offset += snprintf(buffer + offset, buffer_size - offset, "C Flags: %s\n", C_FLAGS);
    offset += snprintf(buffer + offset, buffer_size - offset, "C Release Flags: %s\n", C_FLAGS_RELEASE);
    offset += snprintf(buffer + offset, buffer_size - offset, "C Debug Flags: %s\n", C_FLAGS_DEBUG);
    offset += snprintf(buffer + offset, buffer_size - offset, "Build Host: %s\n", BUILD_HOST);
    offset += snprintf(buffer + offset, buffer_size - offset, "CMake Timestamp: %s\n", CMAKE_TIMESTAMP);
              snprintf(buffer + offset, buffer_size - offset, "System: %s (%s)\n", SYSTEM_NAME, SYSTEM_PROCESSOR);
}

#ifdef __cplusplus
}
#endif

#endif //N_BODY_SIMULATION_BUILD_UTILS_H