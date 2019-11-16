project( cvMat_pybind )
cmake_minimum_required(VERSION 3.15)
set (CMAKE_CXX_STANDARD 11)
set(PYBIND11_CPP_STANDARD -std=c++11)

# specify here which python interpreter you want to compile this module for
set(PYTHON_EXECUTABLE "~/miniconda3/bin/python3")

add_compile_definitions(BUILDING_PYTHON_MODULE) # disable rcpp specific import

# The following requires pybind installed in the compiling system
add_subdirectory(pybind11)
pybind11_add_module(segmentation_algos src/omega.cpp  src/listPoint.cpp src/Point.cpp src/Costs.cpp src/pybind_def.cpp)

add_executable(out src/test.cpp src/omega.cpp src/Point.cpp src/ListPoint.cpp src/Costs.cpp)
