cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

foreach(required_variable
        CIS2M_BUILD_DIR
        CIS2M_SOURCE_DIR
        CIS2M_INSTALL_CMAKEDIR
        CIS2M_CTEST_COMMAND)
    if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
        message(FATAL_ERROR "Required variable ${required_variable} is not set")
    endif()
endforeach()

set(package_test_root "${CIS2M_BUILD_DIR}/package-test")
set(install_prefix "${package_test_root}/install")
set(consumer_build_dir "${package_test_root}/build")

file(REMOVE_RECURSE "${package_test_root}")

set(install_command
    "${CMAKE_COMMAND}"
    --install "${CIS2M_BUILD_DIR}"
    --prefix "${install_prefix}")
if(DEFINED CIS2M_TEST_CONFIG AND NOT "${CIS2M_TEST_CONFIG}" STREQUAL "")
    list(APPEND install_command --config "${CIS2M_TEST_CONFIG}")
endif()

execute_process(
    COMMAND ${install_command}
    RESULT_VARIABLE install_result)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "Installing cis2m for the package test failed")
endif()

set(configure_command
    "${CMAKE_COMMAND}"
    -S "${CIS2M_SOURCE_DIR}/test/cmake"
    -B "${consumer_build_dir}"
    "-Dcis2m_DIR=${install_prefix}/${CIS2M_INSTALL_CMAKEDIR}")
if(DEFINED CIS2M_EIGEN3_DIR AND NOT "${CIS2M_EIGEN3_DIR}" STREQUAL "")
    list(APPEND configure_command "-DEigen3_DIR=${CIS2M_EIGEN3_DIR}")
endif()

execute_process(
    COMMAND ${configure_command}
    RESULT_VARIABLE configure_result)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "Configuring the cis2m package consumer failed")
endif()

set(build_command "${CMAKE_COMMAND}" --build "${consumer_build_dir}")
if(DEFINED CIS2M_TEST_CONFIG AND NOT "${CIS2M_TEST_CONFIG}" STREQUAL "")
    list(APPEND build_command --config "${CIS2M_TEST_CONFIG}")
endif()

execute_process(
    COMMAND ${build_command}
    RESULT_VARIABLE build_result)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "Building the cis2m package consumer failed")
endif()

set(test_command
    "${CMAKE_COMMAND}"
    -E chdir "${consumer_build_dir}"
    "${CIS2M_CTEST_COMMAND}"
    --output-on-failure)
if(DEFINED CIS2M_TEST_CONFIG AND NOT "${CIS2M_TEST_CONFIG}" STREQUAL "")
    list(APPEND test_command -C "${CIS2M_TEST_CONFIG}")
endif()

execute_process(
    COMMAND ${test_command}
    RESULT_VARIABLE test_result)
if(NOT test_result EQUAL 0)
    message(FATAL_ERROR "Running the cis2m package consumer failed")
endif()
