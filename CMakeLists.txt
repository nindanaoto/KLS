cmake_minimum_required(VERSION 3.18)

project(KLS
  VERSION 0.1.0
  DESCRIPTION "Standalone sparse direct solver for SPICE-style repeated linear solves"
  LANGUAGES C)

option(KLS_NATIVE_ARCH "Compile for the build machine's CPU (-march=native)" ON)
if(KLS_NATIVE_ARCH)
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-march=native" KLS_HAVE_MARCH_NATIVE)
  if(KLS_HAVE_MARCH_NATIVE)
    add_compile_options(-march=native)
  endif()
endif()

option(KLS_BUILD_TESTS "Build KLS tests" ON)
option(KLS_BUILD_BENCHMARKS "Build KLS benchmark tools" ON)
option(KLS_BUILD_CKTSO_COMPARE "Build optional CKTSO comparison tool" OFF)
option(KLS_BUILD_KLU2_COMPARE "Build optional Trilinos KLU2 comparison tool" OFF)
option(KLS_ENABLE_METIS "Enable METIS nested-dissection ordering support" ON)
option(KLS_USE_SYSTEM_METIS "Use a system METIS library instead of the bundled submodule" OFF)
option(KLS_ENABLE_SCOTCH "Enable SCOTCH nested-dissection ordering support" ON)
option(KLS_USE_SYSTEM_SCOTCH "Use a system SCOTCH library instead of the bundled submodule" OFF)
option(KLS_ENABLE_SPRAL_SCALING "Enable SPRAL matching/scaling support" ON)
option(KLS_USE_SYSTEM_SPRAL "Use a system SPRAL library instead of the bundled submodule" OFF)
option(KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE "Acknowledge that the selected system SPRAL library is redistributable with LGPL KLS" OFF)
option(KLS_ENABLE_CBLAS_SUPERNODE "Build optional CBLAS kernels for compact row-supernode experiments" OFF)
set(CKTSO_ROOT "" CACHE PATH "Path to a local CKTSO distribution")
set(KLU2_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../Trilinos/packages/amesos2/src/KLU2" CACHE PATH "Path to a Trilinos Amesos2 KLU2 source directory")
set(TEUCHOS_CORE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../Trilinos/packages/rol/src/compatibility/teuchos-lite" CACHE PATH "Path to Trilinos Teuchos-compatible headers for KLU2")
set(KLS_SYSTEM_METIS_INCLUDE_DIR "" CACHE PATH "Path to a system METIS include directory")
set(KLS_SYSTEM_METIS_LIBRARY "" CACHE FILEPATH "Path to a system METIS library")
set(KLS_SYSTEM_SCOTCH_INCLUDE_DIR "" CACHE PATH "Path to a system SCOTCH include directory")
set(KLS_SYSTEM_SCOTCH_LIBRARY "" CACHE FILEPATH "Path to a system SCOTCH library")
set(KLS_SYSTEM_SCOTCHERR_LIBRARY "" CACHE FILEPATH "Path to an optional system SCOTCH error library")
set(KLS_SYSTEM_SPRAL_INCLUDE_DIR "" CACHE PATH "Path to a system SPRAL include directory")
set(KLS_SYSTEM_SPRAL_LIBRARY "" CACHE FILEPATH "Path to a system SPRAL library")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

find_library(KLS_MATH_LIBRARY m)
find_package(Threads REQUIRED)
find_package(Python3 COMPONENTS Interpreter QUIET)
include(CheckIncludeFile)

if(KLS_ENABLE_CBLAS_SUPERNODE)
  find_package(BLAS QUIET)
  check_include_file(cblas.h KLS_HAVE_CBLAS_HEADER)
  if(BLAS_FOUND AND KLS_HAVE_CBLAS_HEADER)
    add_library(kls_cblas INTERFACE)
    target_compile_definitions(kls_cblas
      INTERFACE
        KLS_HAVE_CBLAS=1)
    if(TARGET BLAS::BLAS)
      target_link_libraries(kls_cblas
        INTERFACE
          BLAS::BLAS)
    else()
      target_link_libraries(kls_cblas
        INTERFACE
          ${BLAS_LIBRARIES})
    endif()
  else()
    message(STATUS
      "CBLAS compact row-supernode kernels disabled: cblas.h or BLAS library not found")
  endif()
endif()

if(KLS_ENABLE_METIS)
  if(KLS_USE_SYSTEM_METIS)
    find_path(KLS_SYSTEM_METIS_INCLUDE_DIR
      NAMES metis.h
      PATHS ${KLS_SYSTEM_METIS_INCLUDE_DIR})
    find_library(KLS_SYSTEM_METIS_LIBRARY
      NAMES metis
      PATHS ${KLS_SYSTEM_METIS_LIBRARY})
    if(NOT KLS_SYSTEM_METIS_INCLUDE_DIR OR NOT KLS_SYSTEM_METIS_LIBRARY)
      message(FATAL_ERROR
        "KLS_USE_SYSTEM_METIS=ON requires KLS_SYSTEM_METIS_INCLUDE_DIR and "
        "KLS_SYSTEM_METIS_LIBRARY, or a discoverable system METIS install")
    endif()

    add_library(kls_metis INTERFACE)
    target_include_directories(kls_metis
      INTERFACE
        "${KLS_SYSTEM_METIS_INCLUDE_DIR}")
    target_link_libraries(kls_metis
      INTERFACE
        "${KLS_SYSTEM_METIS_LIBRARY}")
    target_compile_definitions(kls_metis
      INTERFACE
        KLS_HAVE_METIS=1
        IDXTYPEWIDTH=64
        REALTYPEWIDTH=64)
  else()
    set(KLS_METIS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/metis")
    set(KLS_GKLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/gklib")
    if(NOT EXISTS "${KLS_METIS_ROOT}/include/metis.h" OR
       NOT EXISTS "${KLS_GKLIB_ROOT}/include/GKlib.h")
      message(FATAL_ERROR
        "Bundled METIS support requires initialized submodules. Run: "
        "git submodule update --init --recursive third_party/metis third_party/gklib")
    endif()

    file(GLOB KLS_GKLIB_SOURCES CONFIGURE_DEPENDS
      "${KLS_GKLIB_ROOT}/src/*.c")
    add_library(kls_gklib STATIC
      ${KLS_GKLIB_SOURCES})
    set_target_properties(kls_gklib PROPERTIES
      POSITION_INDEPENDENT_CODE ON)
    target_include_directories(kls_gklib
      PUBLIC
        "${KLS_GKLIB_ROOT}/include")
    target_compile_definitions(kls_gklib
      PUBLIC
        LINUX
        _POSIX_C_SOURCE=200809L
        NDEBUG
        NDEBUG2)
    if(KLS_MATH_LIBRARY)
      target_link_libraries(kls_gklib
        PUBLIC
          "${KLS_MATH_LIBRARY}")
    endif()
    target_compile_options(kls_gklib PRIVATE
      $<$<C_COMPILER_ID:GNU,Clang>:-Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-unused-label -Wno-pedantic>)

    file(GLOB KLS_METIS_SOURCES CONFIGURE_DEPENDS
      "${KLS_METIS_ROOT}/libmetis/*.c")
    add_library(kls_metis STATIC
      ${KLS_METIS_SOURCES})
    set_target_properties(kls_metis PROPERTIES
      POSITION_INDEPENDENT_CODE ON)
    target_include_directories(kls_metis
      PUBLIC
        "${KLS_METIS_ROOT}/include"
      PRIVATE
        "${KLS_METIS_ROOT}/libmetis")
    target_compile_definitions(kls_metis
      PUBLIC
        KLS_HAVE_METIS=1
        IDXTYPEWIDTH=64
        REALTYPEWIDTH=64
      PRIVATE
        LINUX
        _FILE_OFFSET_BITS=64)
    target_link_libraries(kls_metis
      PRIVATE
        kls_gklib)
    if(KLS_MATH_LIBRARY)
      target_link_libraries(kls_metis
        PRIVATE
          "${KLS_MATH_LIBRARY}")
    endif()
    target_compile_options(kls_metis PRIVATE
      $<$<C_COMPILER_ID:GNU,Clang>:-Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-unused-label -Wno-pedantic>)
  endif()
endif()

if(KLS_ENABLE_SCOTCH)
  if(KLS_USE_SYSTEM_SCOTCH)
    find_path(KLS_SYSTEM_SCOTCH_INCLUDE_DIR
      NAMES scotch.h
      PATHS ${KLS_SYSTEM_SCOTCH_INCLUDE_DIR})
    find_library(KLS_SYSTEM_SCOTCH_LIBRARY
      NAMES scotch
      PATHS ${KLS_SYSTEM_SCOTCH_LIBRARY})
    find_library(KLS_SYSTEM_SCOTCHERR_LIBRARY
      NAMES scotcherr
      PATHS ${KLS_SYSTEM_SCOTCHERR_LIBRARY})
    if(NOT KLS_SYSTEM_SCOTCH_INCLUDE_DIR OR NOT KLS_SYSTEM_SCOTCH_LIBRARY)
      message(FATAL_ERROR
        "KLS_USE_SYSTEM_SCOTCH=ON requires KLS_SYSTEM_SCOTCH_INCLUDE_DIR and "
        "KLS_SYSTEM_SCOTCH_LIBRARY, or a discoverable system SCOTCH install")
    endif()

    add_library(kls_scotch INTERFACE)
    target_include_directories(kls_scotch
      INTERFACE
        "${KLS_SYSTEM_SCOTCH_INCLUDE_DIR}")
    target_link_libraries(kls_scotch
      INTERFACE
        "${KLS_SYSTEM_SCOTCH_LIBRARY}")
    if(KLS_SYSTEM_SCOTCHERR_LIBRARY)
      target_link_libraries(kls_scotch
        INTERFACE
          "${KLS_SYSTEM_SCOTCHERR_LIBRARY}")
    endif()
    if(KLS_MATH_LIBRARY)
      target_link_libraries(kls_scotch
        INTERFACE
          "${KLS_MATH_LIBRARY}")
    endif()
    target_compile_definitions(kls_scotch
      INTERFACE
        KLS_HAVE_SCOTCH=1)
  else()
    set(KLS_SCOTCH_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/scotch")
    if(NOT EXISTS "${KLS_SCOTCH_ROOT}/CMakeLists.txt")
      message(FATAL_ERROR
        "Bundled SCOTCH support requires initialized submodules. Run: "
        "git submodule update --init --recursive third_party/scotch")
    endif()

    set(BUILD_FORTRAN OFF CACHE BOOL "" FORCE)
    set(BUILD_PTSCOTCH OFF CACHE BOOL "" FORCE)
    set(BUILD_LIBESMUMPS OFF CACHE BOOL "" FORCE)
    set(BUILD_LIBSCOTCHMETIS OFF CACHE BOOL "" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(USE_ZLIB OFF CACHE BOOL "" FORCE)
    set(USE_BZ2 OFF CACHE BOOL "" FORCE)
    set(USE_LZMA OFF CACHE BOOL "" FORCE)
    set(INTSIZE 64 CACHE STRING "" FORCE)
    set(IDXSIZE 64 CACHE STRING "" FORCE)
    set(THREADS ON CACHE BOOL "" FORCE)
    set(SCOTCH_DETERMINISTIC "FIXED_SEED" CACHE STRING "" FORCE)
    add_subdirectory("${KLS_SCOTCH_ROOT}" third_party/scotch EXCLUDE_FROM_ALL)
    add_library(kls_scotch INTERFACE)
    target_link_libraries(kls_scotch
      INTERFACE
        SCOTCH::scotch
        SCOTCH::scotcherr)
    target_compile_definitions(kls_scotch
      INTERFACE
        KLS_HAVE_SCOTCH=1)
    set_target_properties(scotch scotcherr PROPERTIES
      POSITION_INDEPENDENT_CODE ON)
  endif()
endif()

if(KLS_ENABLE_SPRAL_SCALING)
  if(KLS_USE_SYSTEM_SPRAL)
    if(NOT KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE)
      message(FATAL_ERROR
        "KLS_USE_SYSTEM_SPRAL=ON may only be used with SPRAL or MC64-style "
        "code that is redistributable with LGPL-2.1-or-later KLS. After "
        "verifying the selected library's license, reconfigure with "
        "-DKLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON.")
    endif()
    find_path(KLS_SYSTEM_SPRAL_INCLUDE_DIR
      NAMES spral_scaling.h
      PATHS ${KLS_SYSTEM_SPRAL_INCLUDE_DIR})
    find_library(KLS_SYSTEM_SPRAL_LIBRARY
      NAMES spral
      PATHS ${KLS_SYSTEM_SPRAL_LIBRARY})
    if(NOT KLS_SYSTEM_SPRAL_INCLUDE_DIR OR NOT KLS_SYSTEM_SPRAL_LIBRARY)
      message(FATAL_ERROR
        "KLS_USE_SYSTEM_SPRAL=ON requires KLS_SYSTEM_SPRAL_INCLUDE_DIR "
        "and KLS_SYSTEM_SPRAL_LIBRARY, or a discoverable system SPRAL install")
    endif()

    add_library(kls_spral_scaling INTERFACE)
    target_include_directories(kls_spral_scaling
      INTERFACE
        "${KLS_SYSTEM_SPRAL_INCLUDE_DIR}")
    target_link_libraries(kls_spral_scaling
      INTERFACE
        "${KLS_SYSTEM_SPRAL_LIBRARY}")
    target_compile_definitions(kls_spral_scaling
      INTERFACE
        KLS_HAVE_SPRAL_SCALING=1)
  else()
    enable_language(Fortran)
    set(KLS_SPRAL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/spral")
    if(NOT EXISTS "${KLS_SPRAL_ROOT}/include/spral_scaling.h" OR
       NOT EXISTS "${KLS_SPRAL_ROOT}/src/scaling.f90" OR
       NOT EXISTS "${KLS_SPRAL_ROOT}/src/matrix_util.f90" OR
       NOT EXISTS "${KLS_SPRAL_ROOT}/interfaces/C/scaling.f90" OR
       NOT EXISTS "${KLS_SPRAL_ROOT}/LICENCE")
      message(FATAL_ERROR
        "Bundled SPRAL scaling support requires an initialized submodule. Run: "
        "git submodule update --init --recursive third_party/spral")
    endif()
    file(READ "${KLS_SPRAL_ROOT}/LICENCE" KLS_SPRAL_LICENSE_TEXT LIMIT 4096)
    string(FIND "${KLS_SPRAL_LICENSE_TEXT}"
      "Redistribution and use in source and binary forms" KLS_SPRAL_BSD_TEXT)
    if(KLS_SPRAL_BSD_TEXT EQUAL -1)
      message(FATAL_ERROR
        "Bundled SPRAL scaling support did not expose the expected "
        "redistribution-compatible BSD license text in third_party/spral/LICENCE")
    endif()

    add_library(kls_spral_scaling STATIC
      "${KLS_SPRAL_ROOT}/src/matrix_util.f90"
      "${KLS_SPRAL_ROOT}/src/scaling.f90"
      "${KLS_SPRAL_ROOT}/interfaces/C/scaling.f90")
    set_target_properties(kls_spral_scaling PROPERTIES
      POSITION_INDEPENDENT_CODE ON
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/spral_modules")
    target_include_directories(kls_spral_scaling
      INTERFACE
        "${KLS_SPRAL_ROOT}/include")
    target_compile_definitions(kls_spral_scaling
      INTERFACE
        KLS_HAVE_SPRAL_SCALING=1)
  endif()
endif()

set(KLS_SUITESPARSE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/suitesparse")
set(KLS_SUITESPARSE_INCLUDE_DIRS
  "${KLS_SUITESPARSE_ROOT}/KLU/Include"
  "${KLS_SUITESPARSE_ROOT}/AMD/Include"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Include"
  "${KLS_SUITESPARSE_ROOT}/COLAMD/Include"
  "${KLS_SUITESPARSE_ROOT}/CCOLAMD/Include"
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
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_sort.c"
  "${KLS_SUITESPARSE_ROOT}/KLU/Source/trilinos_klu_tsolve.c")
set_source_files_properties(${KLS_KLU_DLONG_SOURCES}
  PROPERTIES COMPILE_DEFINITIONS DLONG)

set(KLS_CAMD_DLONG_SOURCES
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_1.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_aat.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_control.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_defaults.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_dump.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_global.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_info.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_order.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_preprocess.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_valid.c")
set_source_files_properties(${KLS_CAMD_DLONG_SOURCES}
  PROPERTIES COMPILE_DEFINITIONS DLONG)

set(KLS_VENDOR_SOURCES
  ${KLS_KLU_LONG_SOURCES}
  ${KLS_KLU_DLONG_SOURCES}
  ${KLS_CAMD_DLONG_SOURCES}
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
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_l2.c"
  "${KLS_SUITESPARSE_ROOT}/CAMD/Source/trilinos_camd_l_postorder.c"
  "${KLS_SUITESPARSE_ROOT}/CCOLAMD/Source/trilinos_ccolamd_global.c"
  "${KLS_SUITESPARSE_ROOT}/CCOLAMD/Source/trilinos_ccolamd_l.c"
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
    kls_suitesparse
    Threads::Threads)
if(KLS_ENABLE_METIS)
  target_link_libraries(kls
    PRIVATE
      kls_metis)
  # optional prebuilt mt-metis (threaded NodeND experiment, env-gated
  # at runtime via KLS_MT_ND); build it once with its own configure/make
  set(KLS_MTMETIS_LIB
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/mtmetis/build/Linux-x86_64/lib/libmtmetis_hidden.o")
  if(EXISTS "${KLS_MTMETIS_LIB}")
    add_library(kls_mtmetis INTERFACE)
    target_include_directories(kls_mtmetis
      INTERFACE
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/mtmetis/include")
    find_package(OpenMP REQUIRED COMPONENTS C)
    target_link_libraries(kls_mtmetis
      INTERFACE
        "${KLS_MTMETIS_LIB}"
        OpenMP::OpenMP_C)
    target_compile_definitions(kls_mtmetis
      INTERFACE
        KLS_HAVE_MTMETIS=1)
    target_link_libraries(kls
      PRIVATE
        kls_mtmetis)
  endif()
endif()
if(KLS_ENABLE_SCOTCH)
  target_link_libraries(kls
    PRIVATE
      kls_scotch)
endif()
if(KLS_ENABLE_SPRAL_SCALING)
  target_link_libraries(kls
    PRIVATE
      kls_spral_scaling)
endif()
if(TARGET kls_cblas)
  target_link_libraries(kls
    PRIVATE
      kls_cblas)
endif()

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

  find_path(KLS_SYSTEM_KLU_INCLUDE_DIR
    NAMES klu.h
    PATH_SUFFIXES suitesparse)
  find_library(KLS_SYSTEM_KLU_LIBRARY
    NAMES klu)
  if(KLS_SYSTEM_KLU_INCLUDE_DIR AND KLS_SYSTEM_KLU_LIBRARY)
    add_executable(klu_width_compare
      bench/klu_width_compare.c)
    target_include_directories(klu_width_compare
      PRIVATE
        "${KLS_SYSTEM_KLU_INCLUDE_DIR}")
    target_link_libraries(klu_width_compare
      PRIVATE
        "${KLS_SYSTEM_KLU_LIBRARY}")
    if(KLS_MATH_LIBRARY)
      target_link_libraries(klu_width_compare PRIVATE "${KLS_MATH_LIBRARY}")
    endif()
    target_compile_options(klu_width_compare PRIVATE
      $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
  else()
    message(STATUS
      "klu_width_compare disabled: system SuiteSparse KLU headers/library not found")
  endif()
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

if(KLS_BUILD_KLU2_COMPARE)
  if(NOT EXISTS "${KLU2_ROOT}/Include/klu2.h" OR
     NOT EXISTS "${KLU2_ROOT}/Source/klu2_factor.hpp")
    message(FATAL_ERROR
      "KLS_BUILD_KLU2_COMPARE requires -DKLU2_ROOT=/path/to/Trilinos/packages/amesos2/src/KLU2")
  endif()
  if(NOT EXISTS "${TEUCHOS_CORE_ROOT}/Teuchos_ScalarTraits.hpp")
    message(FATAL_ERROR
      "KLS_BUILD_KLU2_COMPARE requires -DTEUCHOS_CORE_ROOT=/path/to/Teuchos-compatible headers")
  endif()
  enable_language(CXX)
  add_executable(klu2_compare
    bench/klu2_compare.cpp)
  target_include_directories(klu2_compare
    PRIVATE
      ${KLS_SUITESPARSE_INCLUDE_DIRS})
  target_include_directories(klu2_compare
    SYSTEM PRIVATE
      "${KLU2_ROOT}/Include"
      "${KLU2_ROOT}/Source"
      "${TEUCHOS_CORE_ROOT}")
  target_link_libraries(klu2_compare PRIVATE kls_suitesparse)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(klu2_compare PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_features(klu2_compare PRIVATE cxx_std_17)
  target_compile_options(klu2_compare PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
endif()

if(KLS_BUILD_TESTS)
  enable_testing()
  if(Python3_Interpreter_FOUND)
    add_test(
      NAME kls_license_boundary
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/scripts/audit_license_boundary.py"
              --root "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
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
