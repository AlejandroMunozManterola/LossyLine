cmake_minimum_required(VERSION 3.12) 

message(STATUS "Creating build system for lossylin...")

add_library(lossylin "lossylin.cpp")

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(lossylin Eigen3::Eigen)

find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(lossylin nlohmann_json nlohmann_json::nlohmann_json)