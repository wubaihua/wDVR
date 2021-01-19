cmake_minimum_required(VERSION 2.8)

project(wDVR)

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER "ifort")
set(CMAKE_Fortran_FLAGS "-mkl")
set(CMAKE_BUILD_TYPE Release)

file(GLOB_RECURSE sources src/*.f90)

add_executable(wDVR ${sources})
