cmake_minimum_required(VERSION 3.12)

message(STATUS "Creating build system for tests...")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(lossylin_tests "TestLossyLin.cpp")
target_link_libraries(lossylin_tests lossylin GTest::gtest GTest::gtest_main)
