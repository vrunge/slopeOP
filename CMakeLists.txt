project( pybind_example )
cmake_minimum_required( VERSION 2.8 )
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -pthread -Wall -Werror")
set(PYBIND11_CPP_STANDARD -std=c++11)
# set(PYTHON_EXECUTABLE "/Users/marco/miniconda3/bin/python3")

add_compile_definitions(BUILDING_PYTHON_MODULE) # disable rcpp specific import

# The following requires pybind installed in the compiling system
find_package(pybind11 REQUIRED)
pybind11_add_module(segmentation_algos src/omega.cpp  src/listPoint.cpp src/Point.cpp src/Costs.cpp src/pybind_def.cpp)


add_executable(out src/test.cpp src/omega.cpp src/Point.cpp src/ListPoint.cpp src/Costs.cpp)
target