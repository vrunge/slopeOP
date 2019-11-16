project( slopeOP_pybind )
cmake_minimum_required(VERSION 3.15)
set (CMAKE_CXX_STANDARD 11)

add_compile_definitions(DISABLE_RCPP) # disable rcpp specific import

# specify here which python interpreter you want to compile this module for
# set(PYTHON_EXECUTABLE "~/miniconda3/bin/python3")

option(BUILD_TESTS "Build tests." ON)
option(BUILD_PYTHON_MODULE "Build python module" OFF)

message("BUILD_TESTS = ${BUILD_TESTS}")

# The following requires pybind installed in the compiling system
if (BUILD_PYTHON_MODULE)
    set(PYBIND11_CPP_STANDARD -std=c++11)
    message(">>> Python module for SlopeOP will be built")
    add_subdirectory(pybind11)
    pybind11_add_module(slopeOP src/omega.cpp  src/listPoint.cpp src/Point.cpp src/Costs.cpp src/pybind_def.cpp)
endif()

if (BUILD_TESTS)
    message(">>> Test will be built")
    add_executable(test src/test.cpp src/omega.cpp src/Point.cpp src/ListPoint.cpp src/Costs.cpp)
endif()