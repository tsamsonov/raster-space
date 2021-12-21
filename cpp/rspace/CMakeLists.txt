cmake_minimum_required(VERSION 2.8.12)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(rspace)

SET (PYBIND11_PYTHON_VERSION 3.9 CACHE STRING "")
SET (CMAKE_CXX_STANDARD 17)
SET (PYTHON_EXECUTABLE "/usr/local/bin/python3.9")
SET (VERSION "")

add_subdirectory(pybind11)
pybind11_add_module(rspace rspace.cpp)
