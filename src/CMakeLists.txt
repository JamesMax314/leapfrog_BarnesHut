cmake_minimum_required(VERSION 3.14)
project(barnesHut)

set(CMAKE_CXX_STANDARD 14)

SET(APP_EXE barnesHut)

#FIND_PACKAGE(pybind11 REQUIRED)

ADD_LIBRARY(trees STATIC trees.cpp)
ADD_LIBRARY(bodies STATIC bodies.cpp)
ADD_LIBRARY(vecMaths STATIC vecMaths.cpp)
ADD_LIBRARY(leapfrog STATIC leapfrog.cpp)
#ADD_LIBRARY(save STATIC save.cpp)
add_executable(barnesHut main.cpp)

include_directories(${PROJECT_SOURCE_DIR})
#include_directories(.//include)
TARGET_LINK_LIBRARIES(${APP_EXE} trees)
TARGET_LINK_LIBRARIES(${APP_EXE} vecMaths)
TARGET_LINK_LIBRARIES(${APP_EXE} leapfrog)
TARGET_LINK_LIBRARIES(${APP_EXE} bodies)
#include_directories(include)
