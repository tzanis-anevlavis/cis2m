cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(cis2m_package_test LANGUAGES CXX)

include(CTest)

find_package(cis2m REQUIRED CONFIG)

add_executable(test_cmake test_cmake.cc)
target_link_libraries(test_cmake PRIVATE cis2m::cis2m)

add_test(NAME test_cmake COMMAND test_cmake)
