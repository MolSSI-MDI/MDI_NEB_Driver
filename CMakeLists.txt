cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)

project(NEB_Driver
        LANGUAGES C CXX)

add_subdirectory(lib)
add_subdirectory(NEB_Driver)
