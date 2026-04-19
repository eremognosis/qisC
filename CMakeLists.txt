cmake_minimum_required(VERSION 4.2)
project(qisC C)

set(CMAKE_C_STANDARD 11)

include_directories(include)
include_directories(include/core)

add_executable(qisC
        include/core/error.h
        include/core/qubit.h
        src/core/qubit.c)
