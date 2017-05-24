cmake_minimum_required(VERSION 3.7)
project(ClangBug)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -O3)
set(SOURCE_FILES main.cpp)
add_executable(ClangBug ${SOURCE_FILES})
