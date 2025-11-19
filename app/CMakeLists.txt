add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/app/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_BINARY_DIR}/generated")

#target_compile_options(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -march=native")

if(USE_MPI)
    find_and_link_mpi_for_target(${PROJECT_NAME} PRIVATE)
endif()

if(USE_OPENMP)
    find_and_link_openmp_for_target(${PROJECT_NAME} PRIVATE)
endif ()

target_link_libraries(${PROJECT_NAME} PRIVATE benchmark)
target_link_libraries(${PROJECT_NAME} PRIVATE utils)
target_link_libraries(${PROJECT_NAME} PRIVATE print_debug)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)



set(MAIN_PRINT_DEBUG_LEVEL ${GLOBAL_PRINT_DEBUG_LEVEL} CACHE STRING "Debug level for main")
define_print_debug_level_for_target(${PROJECT_NAME} ${MAIN_PRINT_DEBUG_LEVEL})
