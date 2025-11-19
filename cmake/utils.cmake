include_guard(DIRECTORY)

add_library(utils INTERFACE)
target_include_directories(utils INTERFACE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(utils INTERFACE "${PROJECT_BINARY_DIR}/generated")

