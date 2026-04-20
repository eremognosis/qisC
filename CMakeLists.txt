cmake_minimum_required(VERSION 3.10)
project(qisC VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include(GNUInstallDirs)

option(QISC_BUILD_SHARED "Build the shared library used by the Python wrapper" ON)
option(QISC_BUILD_STATIC "Build a static library" ON)
option(QISC_BUILD_CLI "Build the qisc command-line executable" ON)

set(QISC_PUBLIC_HEADERS
        include/qisc.h
        include/core/error.h
        include/core/qubit.h
        include/core/gate.h
        include/core/tensor.h
        include/gates/gatelib.h
        include/algos/grover.h)

set(QISC_SOURCES
        src/core/error.c
        src/core/qubit.c
        src/core/gate.c
        src/core/tensor.c
        src/gates/gatelib.c
        src/algos/grover.c)

find_library(QISC_MATH_LIBRARY m)

function(qisc_configure_library target_name)
    target_include_directories(${target_name}
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME qisC
            POSITION_INDEPENDENT_CODE ON)

    if (QISC_MATH_LIBRARY)
        target_link_libraries(${target_name} PUBLIC ${QISC_MATH_LIBRARY})
    endif ()
endfunction()

set(QISC_LIBRARY_TARGETS)
set(QISC_DEFAULT_LINK_TARGET)

if (QISC_BUILD_SHARED)
    add_library(qisC SHARED ${QISC_SOURCES} ${QISC_PUBLIC_HEADERS})
    qisc_configure_library(qisC)
    set_target_properties(qisC PROPERTIES
            VERSION ${PROJECT_VERSION}
            SOVERSION 0)
    list(APPEND QISC_LIBRARY_TARGETS qisC)
    set(QISC_DEFAULT_LINK_TARGET qisC)
endif ()

if (QISC_BUILD_STATIC)
    add_library(qisC_static STATIC ${QISC_SOURCES} ${QISC_PUBLIC_HEADERS})
    qisc_configure_library(qisC_static)
    list(APPEND QISC_LIBRARY_TARGETS qisC_static)
    if (NOT QISC_DEFAULT_LINK_TARGET)
        set(QISC_DEFAULT_LINK_TARGET qisC_static)
    endif ()
endif ()

if (NOT QISC_DEFAULT_LINK_TARGET)
    message(FATAL_ERROR "At least one of QISC_BUILD_SHARED or QISC_BUILD_STATIC must be enabled.")
endif ()

if (QISC_BUILD_CLI)
    add_executable(qisc_cli src/cli/qisc_cli.c)
    target_link_libraries(qisc_cli PRIVATE ${QISC_DEFAULT_LINK_TARGET})
    set_target_properties(qisc_cli PROPERTIES OUTPUT_NAME qisc)
endif ()

install(TARGETS ${QISC_LIBRARY_TARGETS}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if (QISC_BUILD_CLI)
    install(TARGETS qisc_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
