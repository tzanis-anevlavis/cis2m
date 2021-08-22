cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(cis2m 
	LANGUAGES CXX
	VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(ortools REQUIRED)


message(STATUS "Project will be installed to ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

include(GNUInstallDirs)

# Directory layout in the build root
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
	${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
	${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
	${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for headers")
set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

foreach (p LIB BIN INCLUDE CMAKE)
	file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path)
	message(STATUS "Installing ${p} components to ${_path}")
	unset(_path)
endforeach()

set(test_folder "test")
set(src_folder "src")
set(header_folder "include")
set(include_files "include/hpolyhedron.hpp;include/cis_generator.hpp;include/brunovskyform.hpp")

enable_testing()

# Library
add_library(cis2m SHARED)

target_sources(cis2m
	PRIVATE
	${src_folder}/hpolyhedron.cc
	${src_folder}/brunovskyform.cc
	${src_folder}/cis_generator.cc
	)

target_include_directories(cis2m
	PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/${INSTALL_INCLUDEDIR}>
	$<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
	)

target_link_libraries(cis2m
	PUBLIC
	ortools::ortools
	Eigen3::Eigen
	)

target_compile_definitions(cis2m 
	PUBLIC CIS2M_DEBUG
	INTERFACE
	$<INSTALL_INTERFACE:USING_cis2m>
	)

set_target_properties(cis2m
	PROPERTIES
	POSITION_INDEPENDENT_CODE 1
	DEBUG_POSTFIX "_d"
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "cis2m"
	PUBLIC_HEADER "${include_files}" 
	)

install(
	TARGETS cis2m
	EXPORT
	cis2mTargets
	ARCHIVE
	DESTINATION ${INSTALL_LIBDIR}
	COMPONENT lib
	RUNTIME
	DESTINATION ${INSTALL_BINDIR}
	COMPONENT bin
	LIBRARY
	DESTINATION ${INSTALL_DIR}
	COMPONENT lib
	PUBLIC_HEADER
	DESTINATION ${INSTALL_INCLUDEDIR}/cis2m
	COMPONENT dev
	)


# Install rules for the cmake configuration file
install(
	EXPORT
	cis2mTargets
	NAMESPACE
	"cis2m::"
	DESTINATION
	${INSTALL_CMAKEDIR}
	COMPONENT
	dev
	)

include(CMakePackageConfigHelpers)


## Configuration file generation and installation
# 1. Generate file with version information
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/cis2mConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
	)

# 2. Generate cmake config file
configure_package_config_file(
	${PROJECT_SOURCE_DIR}/cmake/cis2mConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/cis2mConfig.cmake
	INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
	)

# 3. Install
install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/cis2mConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cis2mConfigVersion.cmake
	DESTINATION
	${INSTALL_CMAKEDIR}
	)

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
