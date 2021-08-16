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
add_executable(hpoly_affine_test_0 ${test_folder}/test_hpolyhedron_affineT0.cc)
target_link_libraries(hpoly_affine_test_0 cis2m)

add_executable(hpoly_affine_test_1 ${test_folder}/test_hpolyhedron_affineT1.cc)
target_link_libraries(hpoly_affine_test_1 cis2m)

add_executable(hpoly_minkowski_test_0 ${test_folder}/test_hpolyhedron_mink0.cc)
target_link_libraries(hpoly_minkowski_test_0 cis2m)

add_executable(hpoly_minkowski_test_1 ${test_folder}/test_hpolyhedron_mink1.cc)
target_link_libraries(hpoly_minkowski_test_1 cis2m)

add_executable(hpoly_minkowski_test_2 ${test_folder}/test_hpolyhedron_mink2.cc)
target_link_libraries(hpoly_minkowski_test_2 cis2m)

add_executable(hpoly_minkowski_test_3 ${test_folder}/test_hpolyhedron_mink3.cc)
target_link_libraries(hpoly_minkowski_test_3 cis2m)

add_executable(cis_test_0 ${test_folder}/test_cisgenerator0.cc)
target_link_libraries(cis_test_0 cis2m)

add_executable(cis_test_1 ${test_folder}/test_cisgenerator1.cc)
target_link_libraries(cis_test_1 cis2m)


enable_testing()

add_test(
	NAME hpoly_affine_test_0
	COMMAND $<TARGET_FILE:hpoly_affine_test_0>
	)
add_test(
	NAME hpoly_affine_test_1
	COMMAND $<TARGET_FILE:hpoly_affine_test_1>
	)
add_test(
	NAME hpoly_minkowski_test_0
	COMMAND $<TARGET_FILE:hpoly_minkowski_test_0>
	)
add_test(
	NAME hpoly_minkowski_test_1
	COMMAND $<TARGET_FILE:hpoly_minkowski_test_1>
	)
add_test(
	NAME hpoly_minkowski_test_2
	COMMAND $<TARGET_FILE:hpoly_minkowski_test_2>
	)
add_test(
	NAME hpoly_minkowski_test_3
	COMMAND $<TARGET_FILE:hpoly_minkowski_test_3>
	)
