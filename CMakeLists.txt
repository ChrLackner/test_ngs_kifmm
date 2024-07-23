cmake_minimum_required(VERSION 3.28)

include(FetchContent)

project(test_ngs_kifmm)
include(ngsolve_addon.cmake)

#########################################
# FMM3d example
#########################################
FetchContent_Declare(
  fmm3d
  GIT_REPOSITORY https://github.com/flatironinstitute/FMM3D.git
  GIT_TAG master
  EXCLUDE_FROM_ALL
)

add_custom_target(build_fmm3d ALL
  COMMAND make lib -j 9 WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/fmm3d-src
  # BYPRODUCTS ${FETCHCONTENT_BASE_DIR}/kifmm-src/targets/kifmm_rs.h
  USES_TERMINAL
)

add_executable(test_fmm3d src/test_fmm3d.cpp)
add_dependencies(test_fmm3d build_fmm3d)
target_include_directories(test_fmm3d PRIVATE ${FETCHCONTENT_BASE_DIR}/fmm3d-src/c)
target_link_libraries(test_fmm3d PRIVATE ngsolve ${FETCHCONTENT_BASE_DIR}/fmm3d-src/lib/libfmm3d.so)

#########################################
# KIFMM example
#########################################

FetchContent_Declare(
  kifmm
  GIT_REPOSITORY https://github.com/bempp/kifmm.git
  GIT_TAG "enh/c-abi"
  EXCLUDE_FROM_ALL
)

add_custom_target(build_kifmm ALL
  COMMAND cargo build -r WORKING_DIRECTORY ${FETCHCONTENT_BASE_DIR}/kifmm-src
  BYPRODUCTS ${FETCHCONTENT_BASE_DIR}/kifmm-src/targets/kifmm_rs.h
  USES_TERMINAL
)

add_executable(test_kifmm src/test_kifmm.cpp)
target_include_directories(test_kifmm PRIVATE ${FETCHCONTENT_BASE_DIR}/kifmm-src/target)
target_link_libraries(test_kifmm PRIVATE ngsolve ${FETCHCONTENT_BASE_DIR}/kifmm-src/target/release/libkifmm.so)
add_dependencies(test_kifmm build_kifmm)

