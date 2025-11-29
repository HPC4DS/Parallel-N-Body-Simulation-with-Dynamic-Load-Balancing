file(GLOB SANDBOX_SOURCES *.cpp *.c)

set(SANDBOX_PRINT_DEBUG_LEVEL ${GLOBAL_PRINT_DEBUG_LEVEL} CACHE STRING "Debug level for main")

foreach(src ${SANDBOX_SOURCES})
    get_filename_component(exec ${src} NAME_WE)

    add_executable(${exec} ${src})
    find_and_link_mpi_for_target(${exec} PRIVATE)
    target_link_libraries(${exec} PRIVATE benchmark)
    target_link_libraries(${exec} PRIVATE utils)
    target_link_libraries(${exec} PRIVATE print_debug)


    # optional: extra flags
    # target_compile_definitions(${exec} PRIVATE EXPERIMENTAL_MODE)
target_compile_options(${exec} PRIVATE -O3 -ftree-vectorize -march=native)
    target_include_directories(${exec} PUBLIC "${PROJECT_BINARY_DIR}/generated")
    target_include_directories(${exec} PUBLIC ${CMAKE_SOURCE_DIR}/include)

    define_print_debug_level_for_target(${exec} ${SANDBOX_PRINT_DEBUG_LEVEL})
endforeach()