cmake_minimum_required(VERSION 3.7)
project(ClangBug)

set(CMAKE_CXX_STANDARD 11)
#target_compile_options(ClangBug1 PRIVATE -O1)
set(SOURCE_FILES main.cpp)
add_executable(ClangBug0 ${SOURCE_FILES})
add_executable(ClangBug1 ${SOURCE_FILES})
add_executable(ClangBug2 ${SOURCE_FILES})
add_executable(ClangBug3 ${SOURCE_FILES})
set_target_properties( ClangBug3 PROPERTIES COMPILE_FLAGS "-O3" )
set_target_properties( ClangBug1 PROPERTIES COMPILE_FLAGS "-O1" )
set_target_properties( ClangBug2 PROPERTIES COMPILE_FLAGS "-O2" )
set_target_properties( ClangBug0 PROPERTIES COMPILE_FLAGS "-O0" )


