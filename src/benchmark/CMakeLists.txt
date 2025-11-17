add_library(benchmark benchmark.c)

include(${CMAKE_SOURCE_DIR}/cmake/PrintDebug.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/ExternalDependencies.cmake)


# Public headers
target_include_directories(benchmark PUBLIC ${CMAKE_SOURCE_DIR}/include/)

find_and_link_mpi_for_target(benchmark)
target_link_libraries(benchmark PUBLIC print_debug)