cmake_minimum_required(VERSION 2.8.12)
project(rspace)

SET (PYBIND11_PYTHON_VERSION 3.7 CACHE STRING "")
SET (PYTHON_EXECUTABLE "/usr/local/Cellar/python/3.7.2_1/bin/python3.7")
add_subdirectory(pybind11)
pybind11_add_module(rspace rspace.cpp)
