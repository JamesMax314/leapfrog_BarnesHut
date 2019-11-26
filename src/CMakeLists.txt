cmake_minimum_required(VERSION 3.14)
project(barnesHut)

set(CMAKE_CXX_STANDARD 14)

add_executable(barnesHut main.cpp leapfrog.cpp trees.cpp vecMaths.cpp)

include_directories(${PROJECT_SOURCE_DIR})
#include_directories(include)
