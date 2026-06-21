cmake_minimum_required(VERSION 3.18)

project(KLS
  VERSION 0.1.0
  DESCRIPTION "Standalone sparse direct solver for SPICE-style repeated linear solves"
  LANGUAGES C)

option(KLS_BUILD_TESTS "Build KLS tests" ON)
option(KLS_BUILD_BENCHMARKS "Build KLS benchmark tools" ON)
option(KLS_BUILD_CKTSO_COMPARE "Build optional CKTSO comparison tool" OFF)
set(CKTSO_ROOT "" CACHE PATH "Path to a local CKTSO distribution")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

find_library(KLS_MATH_LIBRARY m)

set(KLS_SUITESPARSE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/suitesparse")
set(KLS_SUITESPARSE_INCLUDE_DIRS
  "${KLS_SUITESPARSE_ROOT}/KLU/Include"
  "${KLS_SUITESPARSE_ROOT}/AMD/Include"
  "${KLS_SUITESPARSE_ROOT}/COLAMD/Include"
  "${KLS_SUITESPARSE_ROOT}/BTF/Include"
  "${KLS_SUITESPARSE_ROOT}/UFconfig")

set(KLS_KLU_LONG_SOURCES
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_analyze.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_analyze_given.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_defaults.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_factor.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_free_numeric.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_free_symbolic.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_kernel.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_memory.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_scale.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_l_solve.c")

set(KLS_KLU_DLONG_SOURCES
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_diagnostics.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_refactor.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_tsolve.c")
set_source_files_properties(${KLS_KLU_DLONG_SOURCES}
  PROPERTIES COMPILE_DEFINITIONS DLONG)

set(KLS_VENDOR_SOURCES
  ${KLS_KLU_LONG_SOURCES}
  ${KLS_KLU_DLONG_SOURCES}
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_1.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_2.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_aat.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_control.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_defaults.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_dump.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_global.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_info.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l1.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l2.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_aat.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_order.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_post_tree.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_postorder.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_preprocess.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_l_valid.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_order.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_post_tree.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_postorder.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_preprocess.c"
  "${KLS_SUITESPARSE_ROOT}/AMD/Source/trilinos_amd_valid.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_l_maxtrans.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_l_order.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_l_strongcomp.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_maxtrans.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_order.c"
  "${KLS_SUITESPARSE_ROOT}/BTF/Source/trilinos_btf_strongcomp.c"
  "${KLS_SUITESPARSE_ROOT}/COLAMD/Source/trilinos_colamd.c"
  "${KLS_SUITESPARSE_ROOT}/COLAMD/Source/trilinos_colamd_global.c"
  "${KLS_SUITESPARSE_ROOT}/COLAMD/Source/trilinos_colamd_l.c")

add_library(kls_suitesparse STATIC
  ${KLS_VENDOR_SOURCES})
set_target_properties(kls_suitesparse PROPERTIES
  POSITION_INDEPENDENT_CODE ON)
target_include_directories(kls_suitesparse
  PUBLIC
    ${KLS_SUITESPARSE_INCLUDE_DIRS})
if(KLS_MATH_LIBRARY)
  target_link_libraries(kls_suitesparse
    PUBLIC
      "${KLS_MATH_LIBRARY}")
endif()

add_library(kls
  src/kls.c)

add_library(KLS::kls ALIAS kls)

target_include_directories(kls
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(kls
  PRIVATE
    kls_suitesparse)

target_compile_options(kls PRIVATE
  $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)

if(KLS_BUILD_BENCHMARKS)
  add_executable(kls_bench
    bench/kls_bench.c)
  target_link_libraries(kls_bench PRIVATE KLS::kls)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(kls_bench PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(kls_bench PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
endif()

if(KLS_BUILD_CKTSO_COMPARE)
  if(NOT CKTSO_ROOT)
    message(FATAL_ERROR "KLS_BUILD_CKTSO_COMPARE requires -DCKTSO_ROOT=/path/to/cktso")
  endif()
  enable_language(CXX)
  find_library(CKTSO_LIBRARY
    NAMES cktso
    PATHS
      "${CKTSO_ROOT}/rocky8_x64_gcc850"
      "${CKTSO_ROOT}/centos7_x64_gcc485"
      "${CKTSO_ROOT}/centos6_x64_gcc482"
    NO_DEFAULT_PATH)
  if(NOT CKTSO_LIBRARY)
    message(FATAL_ERROR "Could not find libcktso under CKTSO_ROOT")
  endif()
  add_executable(cktso_compare
    bench/cktso_compare.cpp)
  target_include_directories(cktso_compare PRIVATE "${CKTSO_ROOT}/include")
  target_link_libraries(cktso_compare PRIVATE "${CKTSO_LIBRARY}" m)
  target_compile_features(cktso_compare PRIVATE cxx_std_17)
  target_compile_options(cktso_compare PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
endif()

if(KLS_BUILD_TESTS)
  enable_testing()
  add_executable(kls_smoke
    tests/kls_smoke.c)
  target_link_libraries(kls_smoke PRIVATE KLS::kls)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(kls_smoke PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(kls_smoke PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
  add_test(NAME kls_smoke COMMAND kls_smoke)
endif()

install(TARGETS kls
  EXPORT KLSTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

install(DIRECTORY include/ DESTINATION include)
