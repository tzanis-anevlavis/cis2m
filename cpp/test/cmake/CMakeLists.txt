cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(test_cmake
	LANGUAGES CXX
	VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#find_package(Eigen3 REQUIRED NO_MODULE)
#find_package(ortools REQUIRED)

find_package(cis2m REQUIRED NO_MODULE)
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(test_cmake 
	test_cmake.cc
	)

target_link_libraries(test_cmake 
	Eigen3::Eigen
	cis2m::cis2m
	)
