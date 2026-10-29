cmake_minimum_required(VERSION 3.16)
project(mma_cpp VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

add_library(mma_core
    mma.cpp
    toy2.cpp
    linear_solver.cpp
)
target_include_directories(mma_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(mma_core PUBLIC cxx_std_17)
target_compile_options(mma_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(mma_demo main.cpp)
target_link_libraries(mma_demo PRIVATE mma_core)

if(BUILD_TESTING)
    add_executable(linear_solver_test tests/linear_solver_test.cpp)
    target_link_libraries(linear_solver_test PRIVATE mma_core)
    add_test(NAME linear_solver COMMAND linear_solver_test)
    add_test(NAME toy_problem COMMAND mma_demo)
    set_tests_properties(toy_problem PROPERTIES
        PASS_REGULAR_EXPRESSION "KKT norm:"
    )
endif()
