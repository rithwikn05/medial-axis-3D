cmake_minimum_required(VERSION 3.16)
project(medial_axis_3d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(inspect src/app/inspect.cpp)
target_include_directories(inspect PRIVATE src)

add_executable(test_runner tests/test_runner.cpp)
target_include_directories(test_runner PRIVATE src)
