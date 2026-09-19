cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(cis2m
    VERSION 1.0.0
    LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(CTest)
include(GNUInstallDirs)

find_package(Eigen3 REQUIRED CONFIG)
find_package(ortools REQUIRED CONFIG)

message(STATUS "Project will be installed to ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

set(
    CIS2M_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    CACHE PATH
    "Installation directory for cis2m CMake package files")

set(CIS2M_PUBLIC_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/hpolyhedron.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/cis_generator.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/brunovskytransformation.hpp")

function(cis2m_set_output_directories target_name)
    set_target_properties(
        ${target_name}
        PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
            LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
            RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
endfunction()

add_library(cis2m SHARED)
add_library(cis2m::cis2m ALIAS cis2m)

target_sources(cis2m
    PRIVATE
        src/hpolyhedron.cc
        src/brunovskytransformation.cc
        src/cis_generator.cc
        ${CIS2M_PUBLIC_HEADERS})

target_compile_features(cis2m PUBLIC cxx_std_14)

target_include_directories(cis2m
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(cis2m
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        ortools::ortools)

# Eigen 5 + AppleClang 12 can fail on unavailable ARM complex NEON intrinsics
# such as vcmla*. Disable Eigen vectorization for affected compiler versions.
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
    message(WARNING
        "Detected AppleClang ${CMAKE_CXX_COMPILER_VERSION}. "
        "Enabling EIGEN_DONT_VECTORIZE compatibility fallback.")
    target_compile_definitions(cis2m PUBLIC EIGEN_DONT_VECTORIZE)
endif()

set_target_properties(cis2m
    PROPERTIES
        DEBUG_POSTFIX "_d"
        OUTPUT_NAME "cis2m"
        PUBLIC_HEADER "${CIS2M_PUBLIC_HEADERS}"
        SOVERSION "${PROJECT_VERSION_MAJOR}"
        VERSION "${PROJECT_VERSION}")

cis2m_set_output_directories(cis2m)

install(
    TARGETS cis2m
    EXPORT cis2mTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cis2m" COMPONENT dev)

install(
    EXPORT cis2mTargets
    NAMESPACE cis2m::
    DESTINATION "${CIS2M_INSTALL_CMAKEDIR}"
    COMPONENT dev)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cis2mConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/cis2mConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cis2mConfig.cmake"
    INSTALL_DESTINATION "${CIS2M_INSTALL_CMAKEDIR}")

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/cis2mConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cis2mConfigVersion.cmake"
    DESTINATION "${CIS2M_INSTALL_CMAKEDIR}"
    COMPONENT dev)

if(BUILD_TESTING)
    find_package(GTest 1.14 CONFIG QUIET)

    if(NOT GTest_FOUND)
        include(FetchContent)

        if(POLICY CMP0135)
            cmake_policy(SET CMP0135 NEW)
        endif()

        set(INSTALL_GTEST OFF CACHE BOOL "Do not install GoogleTest with cis2m" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "Do not build unused GoogleMock targets" FORCE)
        FetchContent_Declare(
            gtest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
            URL_HASH SHA256=1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4)
        FetchContent_MakeAvailable(gtest)
    endif()

    function(cis2m_add_unit_test target_name source_file)
        add_executable(${target_name} "${source_file}")
        target_link_libraries(${target_name} PRIVATE cis2m GTest::gtest_main)
        cis2m_set_output_directories(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    cis2m_add_unit_test(test_cisgenerator test/test_cisgenerator.cc)
    cis2m_add_unit_test(test_hpolyhedron test/test_hpolyhedron.cc)
    cis2m_add_unit_test(
        test_brunovskytransformation
        test/test_brunovskytransformation.cc)

    add_test(
        NAME test_cmake_package
        COMMAND "${CMAKE_COMMAND}"
            "-DCIS2M_BUILD_DIR=${PROJECT_BINARY_DIR}"
            "-DCIS2M_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
            "-DCIS2M_INSTALL_CMAKEDIR=${CIS2M_INSTALL_CMAKEDIR}"
            "-DCIS2M_EIGEN3_DIR=${Eigen3_DIR}"
            "-DCIS2M_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
            "-DCIS2M_TEST_CONFIG=$<CONFIG>"
            -P "${PROJECT_SOURCE_DIR}/cmake/RunPackageTest.cmake")
    set_tests_properties(test_cmake_package PROPERTIES LABELS package)
endif()
