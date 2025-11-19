add_executable(
        benchmark_tests
        benchmark_tests.cpp

)
target_link_libraries(
        benchmark_tests
        GTest::gtest_main
)


include(GoogleTest)
gtest_discover_tests(benchmark_tests)