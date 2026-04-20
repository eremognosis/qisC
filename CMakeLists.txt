cmake_minimum_required(VERSION 4.2)
project(qisC C)

set(CMAKE_C_STANDARD 11)

include_directories(include)
include_directories(include/core)

add_library(qisC STATIC
        include/core/error.h
        src/core/error.c
        include/core/qubit.h
        src/core/qubit.c
        include/core/gate.h
        src/core/gate.c
        include/core/tensor.h
        src/core/tensor.c
        include/gates/gatelib.h
        src/gates/gatelib.c)

target_link_libraries(qisC PUBLIC m)
