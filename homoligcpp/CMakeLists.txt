cmake_minimum_required(VERSION 3.4)
project(homoligcpp LANGUAGES CXX)

add_subdirectory(pybind11)
pybind11_add_module(homoligcpp aamatrix.cpp aavector.cpp homolig.cpp homolig_bind.cpp)
