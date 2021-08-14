cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(cis2m LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
find_package(ortools REQUIRED)

set(test_folder "test")
set(src_folder "src")

# Library
add_library(cis2m 
	${src_folder}/hpolyhedron.cc
	${src_folder}/brunovskyform.cc
	${src_folder}/cis_generator.cc)

target_include_directories(cis2m
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src
	PUBLIC ${EIGEN3_INCLUDE_DIRS}
	)
target_link_libraries(cis2m
	ortools::ortools
	)

target_compile_definitions(cis2m PUBLIC CIS2M_DEBUG)


# Tests
add_executable(hpoly_test ${test_folder}/test_hpolyhedron.cc)
target_link_libraries(hpoly_test cis2m)

# Tests
add_executable(cis_test ${test_folder}/test_cisgenerator.cc)
target_link_libraries(cis_test cis2m)



enable_testing()

add_test(
	NAME hpoly_test
	COMMAND $<TARGET_FILE:hpoly_test>
	)

