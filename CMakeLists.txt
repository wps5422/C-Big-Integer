cmake_minimum_required(VERSION 3.29)
project(biginteger)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(biginteger integer.hpp
        main.cpp)
