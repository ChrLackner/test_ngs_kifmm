cmake_minimum_required(VERSION 3.28)

include(FetchContent)

project(test_ngs_kifmm)
include(ngsolve_addon.cmake)

FetchContent_Declare(
  kifmm
  GIT_REPOSITORY https://github.com/bempp/kifmm.git
  GIT_TAG "enh/c-abi"
  EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(kifmm)

add_custom_target(build_kifmm ALL
  COMMAND cargo build -r WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/kifmm-src
  BYPRODUCTS ${FETCHCONTENT_BASE_DIR}/kifmm-src/targets/kifmm_rs.h
  USES_TERMINAL
)

add_executable(testfmm src/testfmm.cpp)
target_include_directories(testfmm PRIVATE ${FETCHCONTENT_BASE_DIR}/kifmm-src/target)
target_link_libraries(testfmm PRIVATE ngsolve ${FETCHCONTENT_BASE_DIR}/kifmm-src/target/release/libkifmm.so)
