cmake_minimum_required(VERSION 3.14)
#project(barnesHut)
add_definitions( "-fpic" ) # Linux only
project(treecode)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

#SET(APP_EXE barnesHut)
SET(APP_PY treecode)
set(PYBIND11_PYTHON_VERSION 3.6)
add_subdirectory(.//Include//pybind11)

ADD_LIBRARY(trees STATIC trees.cpp)
ADD_LIBRARY(bodies STATIC bodies.cpp)
ADD_LIBRARY(vecMaths STATIC vecMaths.cpp)
ADD_LIBRARY(leapfrog STATIC leapfrog.cpp)
ADD_LIBRARY(treeShow STATIC treeShow.cpp)
pybind11_add_module(treecode pyInterface.cpp)
#add_executable(APP_EXE main.cpp)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(.//Include)
#include_directories(C://ProgramData//Anaconda3//include) # Needs changing for each machine
#include_directories(/home/james/anaconda3/include/.) # Needs changing for each machine
#TARGET_LINK_LIBRARIES(${APP_EXE} trees vecMaths leapfrog bodies)
TARGET_LINK_LIBRARIES(treecode PRIVATE trees vecMaths leapfrog bodies treeShow)
#include_directories(include)
