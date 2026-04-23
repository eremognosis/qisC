cmake_minimum_required(VERSION 3.10)
project(qisC VERSION 0.1.0 LANGUAGES C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include(GNUInstallDirs)

option(QISC_BUILD_SHARED "Build the shared library" ON)
option(QISC_BUILD_STATIC "Build a static library" ON)
option(QISC_BUILD_CLI    "Build the qisc CLI" ON)

set(QISC_SOURCES
        src/core/error.c
        src/core/qubit.c
        src/core/gate.c
        src/core/tensor.c
        src/gates/gatelib.c
        src/algos/grover.c
        src/algos/fourier.c
)

# Actual header files (The stuff that tells people how to use the stuff)
set(QISC_HEADERS
        include/qisc.h
        include/core/error.h
        include/core/qubit.h
        include/core/gate.h
        include/core/tensor.h
        include/gates/gatelib.h
        include/algos/grover.h
        include/algos/classred.h
        include/algos/fourier.h
        include/algos/phaseest.h
)


find_library(QISC_MATH_LIBRARY m)

# Target configuration helper to stop repeating yourself
function(qisc_configure_target target_name)
    target_include_directories(${target_name}
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    if(QISC_MATH_LIBRARY)
        target_link_libraries(${target_name} PUBLIC ${QISC_MATH_LIBRARY})
    endif()
endfunction()

# --- Library Targets ---

set(QISC_LIB_TARGETS)

if(QISC_BUILD_SHARED)
    add_library(qisC SHARED ${QISC_SOURCES})
    qisc_configure_target(qisC)
    set_target_properties(qisC PROPERTIES
            VERSION ${PROJECT_VERSION}
            SOVERSION 0
            POSITION_INDEPENDENT_CODE ON
    )
    list(APPEND QISC_LIB_TARGETS qisC)
endif()

if(QISC_BUILD_STATIC)
    add_library(qisC_static STATIC ${QISC_SOURCES})
    qisc_configure_target(qisC_static)
    set_target_properties(qisC_static PROPERTIES OUTPUT_NAME qisC)
    list(APPEND QISC_LIB_TARGETS qisC_static)
endif()

# --- CLI Target ---

if(QISC_BUILD_CLI)
    add_executable(qisc_cli src/cli/qisc_cli.c)

    # Link against whatever library we just built
    if(TARGET qisC)
        target_link_libraries(qisc_cli PRIVATE qisC)
    else()
        target_link_libraries(qisc_cli PRIVATE qisC_static)
    endif()

    set_target_properties(qisc_cli PROPERTIES OUTPUT_NAME qisc)
endif()

# --- Installation (For people who actually finish projects) ---

install(TARGETS ${QISC_LIB_TARGETS}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QISC_BUILD_CLI)
    install(TARGETS qisc_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(FILES ${QISC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qisC)
# Also install the directory structure for core/algos/gates
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})