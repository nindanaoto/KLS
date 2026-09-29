cmake_minimum_required(VERSION 3.18)

project(KLS
  VERSION 0.1.0
  DESCRIPTION "Standalone sparse direct solver for SPICE-style repeated linear solves"
  LANGUAGES C)

# Profiles are accepted only by the exact performance implementation that
# produced them. This fingerprint changes with solver sources and codegen.
file(GLOB KLS_TUNING_FINGERPRINT_INPUTS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_accuracy_*.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_componentwise.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_componentwise.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_*.inc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_tuning_profile.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/kls_tuning_profile.h")
list(APPEND KLS_TUNING_FINGERPRINT_INPUTS
  "${CMAKE_CURRENT_SOURCE_DIR}/bench/kls_snb_calibrate.c")
string(TOUPPER "${CMAKE_BUILD_TYPE}" KLS_BUILD_TYPE_UPPER)
set(KLS_TUNING_FINGERPRINT_TEXT
  "${CMAKE_C_COMPILER_ID};${CMAKE_C_COMPILER_VERSION};${CMAKE_C_COMPILER}"
  ";${CMAKE_C_FLAGS};${CMAKE_C_FLAGS_${KLS_BUILD_TYPE_UPPER}};${CMAKE_BUILD_TYPE}"
  ";${CMAKE_SYSTEM_NAME};${CMAKE_SYSTEM_PROCESSOR}")
foreach(KLS_TUNING_INPUT IN LISTS KLS_TUNING_FINGERPRINT_INPUTS)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${KLS_TUNING_INPUT}")
  file(SHA256 "${KLS_TUNING_INPUT}" KLS_TUNING_INPUT_HASH)
  string(APPEND KLS_TUNING_FINGERPRINT_TEXT ";${KLS_TUNING_INPUT_HASH}")
endforeach()
option(KLS_NATIVE_ARCH "Compile for the build machine's CPU (-march=native)" ON)
if(KLS_NATIVE_ARCH)
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-march=native" KLS_HAVE_MARCH_NATIVE)
  if(KLS_HAVE_MARCH_NATIVE)
    add_compile_options(-march=native)
  endif()
endif()

option(KLS_ENABLE_HOSTED_SIMD "Enable compiler SIMD pragmas and FMA contraction in KLS kernels" ON)

option(KLS_BUILD_TESTS "Build KLS tests" ON)
option(KLS_BUILD_BENCHMARKS "Build KLS benchmark tools" ON)
option(KLS_BUILD_AUTOTUNER "Build the host CPU tuning tool" ON)
option(KLS_BUILD_CKTSO_COMPARE "Build optional CKTSO comparison tool" OFF)
option(KLS_BUILD_SUBTREELU_COMPARE "Build optional SubtreeLU comparison tool" OFF)
option(KLS_USE_SYSTEM_METIS "Use a system METIS library instead of the bundled submodule" OFF)
option(KLS_ENABLE_SCOTCH "Enable SCOTCH nested-dissection ordering support" ON)
option(KLS_USE_SYSTEM_SCOTCH "Use a system SCOTCH library instead of the bundled submodule" OFF)
option(KLS_ENABLE_SPRAL_SCALING "Enable SPRAL matching/scaling support" ON)
option(KLS_USE_SYSTEM_SPRAL "Use a system SPRAL library instead of the bundled submodule" OFF)
option(KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE "Acknowledge that the selected system SPRAL library is redistributable with LGPL KLS" OFF)
string(APPEND KLS_TUNING_FINGERPRINT_TEXT
  ";native=${KLS_NATIVE_ARCH};simd=${KLS_ENABLE_HOSTED_SIMD}"
  ";metis=${KLS_USE_SYSTEM_METIS};scotch=${KLS_ENABLE_SCOTCH}"
  ";spral=${KLS_ENABLE_SPRAL_SCALING}")
string(SHA256 KLS_TUNING_BUILD_ID "${KLS_TUNING_FINGERPRINT_TEXT}")
set(CKTSO_ROOT "" CACHE PATH "Path to a local CKTSO distribution")
set(SUBTREELU_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../SubtreeLU" CACHE PATH
  "Path to a local SubtreeLU distribution")
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

# The final caller-system certificate must not silently become a heuristic on
# unsupported compilers or under unsafe reassociation. Compile/link the actual
# checker early, before building the bundled ordering libraries.
include(CheckCSourceCompiles)
unset(KLS_HAVE_ACCURACY_CERTIFICATE CACHE)
set(CMAKE_REQUIRED_LIBRARIES_SAVED "${CMAKE_REQUIRED_LIBRARIES}")
if(KLS_MATH_LIBRARY)
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${KLS_MATH_LIBRARY}")
endif()
check_c_source_compiles("
#include \"${CMAKE_CURRENT_SOURCE_DIR}/src/kls_accuracy_certificate.c\"
int main(void) { double bound;
  return kls_accuracy_certify(0, 0, 0, 0, 0, 0, 0, &bound); }
" KLS_HAVE_ACCURACY_CERTIFICATE)
set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES_SAVED}")
unset(CMAKE_REQUIRED_LIBRARIES_SAVED)
if(NOT KLS_HAVE_ACCURACY_CERTIFICATE)
  message(FATAL_ERROR
    "KLS accuracy verification requires IEEE binary64, compiler-supported IEEE "
    "binary128 (__float128), 128-bit integers, and strict floating-point "
    "semantics. Fast-math and finite-only modes are unsupported. "
    "GCC 14 or later is recommended.")
endif()

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
    # USE_GKRAND: thread-local MT19937 in gklib random.c.  Concurrent
    # METIS calls (raced/speculative analyzes, deterministic parallel
    # NodeNDP leaves) corrupt each other's streams under the stock
    # process-global srand/rand, making orderings run-to-run
    # nondeterministic; with MT19937 each call seeds its own
    # thread-local stream at entry, so every ordering is deterministic
    # under any interleaving.
    option(KLS_GKRAND "thread-local MT19937 RNG for bundled GKlib" ON)
    if(KLS_GKRAND)
      target_compile_definitions(kls_gklib PUBLIC USE_GKRAND)
    endif()
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
    # GCC 15's loop vectorizer changes the resolution of exactly tied
    # Hungarian/auction updates in the bundled SPRAL matcher.  The resulting
    # assignment is valid, but its downstream sparse ordering can retain
    # materially different fill from the same input on different compiler
    # generations.  Keep the rest of SPRAL optimized while preserving the
    # scalar update order used by GCC 13 and earlier.  This is target-local;
    # system SPRAL libraries and every KLS numeric kernel are unaffected.
    target_compile_options(kls_spral_scaling PRIVATE
      $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fno-tree-loop-vectorize>)
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
# Pin hot-loop alignment in the factorization kernels: the serial
# scatter loop's wall time swings up to 15% with link-order code
# layout (identical binaries for the kernel TU measured 13.1 vs 15.1s
# on mac_econ when an unrelated kls.c commit shifted addresses).
target_compile_options(kls_suitesparse PRIVATE
  $<$<C_COMPILER_ID:GNU,Clang>:-falign-functions=64 -falign-loops=32
    -ffp-contract=fast>)

target_include_directories(kls_suitesparse
  PUBLIC
    ${KLS_SUITESPARSE_INCLUDE_DIRS})
if(KLS_MATH_LIBRARY)
  target_link_libraries(kls_suitesparse
    PUBLIC
      "${KLS_MATH_LIBRARY}")
endif()

add_library(kls
  src/kls.c
  src/kls_accuracy_certificate.c
  src/kls_componentwise.c
  src/kls_accuracy_lattice.c
  src/kls_accuracy_recovery.c
  src/kls_accuracy_schedule.c
  src/kls_caller_recovery.inc
  src/kls_lattice_recovery.inc
  src/kls_original_operator.inc
  src/kls_strict_interval.inc
  src/kls_tuning_profile.c
  src/kls_compact_solve.inc
  src/kls_egraph_refactor.inc
  src/kls_first_factor.inc
  src/kls_numeric_engines.inc
  src/kls_numeric_prepare_repair.inc
  src/kls_numeric_schedule.inc
  src/kls_ordering_policy.inc
  src/kls_pts.inc
  src/kls_row_refactor.inc
  src/kls_tuning.inc)

target_compile_definitions(kls PRIVATE
  KLS_TUNING_BUILD_FINGERPRINT="${KLS_TUNING_BUILD_ID}")

# These implementation fragments stay in kls.c's translation unit so moving
# the internal engines does not introduce call boundaries or inhibit inlining.
set_source_files_properties(
  src/kls_caller_recovery.inc
  src/kls_lattice_recovery.inc
  src/kls_original_operator.inc
  src/kls_strict_interval.inc
  src/kls_compact_solve.inc
  src/kls_egraph_refactor.inc
  src/kls_first_factor.inc
  src/kls_numeric_engines.inc
  src/kls_numeric_prepare_repair.inc
  src/kls_numeric_schedule.inc
  src/kls_ordering_policy.inc
  src/kls_pts.inc
  src/kls_row_refactor.inc
  src/kls_tuning.inc
  PROPERTIES HEADER_FILE_ONLY ON)

add_library(KLS::kls ALIAS kls)

target_include_directories(kls
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(kls
  PRIVATE
    kls_suitesparse
    Threads::Threads)
target_link_libraries(kls
  PRIVATE
    kls_metis)
# Optional prebuilt mt-metis supplies threaded separator construction for
# the deterministic NDP forest; build it once with its own configure/make.
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

if(KLS_BUILD_AUTOTUNER)
  add_executable(kls-autotune bench/kls_autotune.c bench/kls_snb_calibrate.c)
  target_include_directories(kls-autotune PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(kls-autotune PRIVATE KLS::kls)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(kls-autotune PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(kls-autotune PRIVATE
    $<TARGET_PROPERTY:kls,COMPILE_OPTIONS>)
  target_compile_definitions(kls-autotune PRIVATE
    $<TARGET_PROPERTY:kls,COMPILE_DEFINITIONS>)
  target_include_directories(kls-autotune PRIVATE
    $<TARGET_PROPERTY:kls,INCLUDE_DIRECTORIES>)
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
target_compile_options(kls PRIVATE
  $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
if(KLS_ENABLE_HOSTED_SIMD)
  target_compile_options(kls PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-fopenmp-simd -ffp-contract=fast>)
endif()

if(KLS_BUILD_BENCHMARKS)
  add_executable(kls_bench
    bench/kls_bench.c)
  target_link_libraries(kls_bench PRIVATE KLS::kls)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(kls_bench PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(kls_bench PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)

  # Always provide a stock-KLU baseline from the exact vendored sources used
  # by KLS.  This avoids making the paper benchmark depend on whichever
  # SuiteSparse version happens to be installed on the host.
  add_executable(klu_compare
    bench/klu_width_compare.c)
  target_compile_definitions(klu_compare PRIVATE KLS_VENDORED_KLU_ONLY=1)
  target_link_libraries(klu_compare PRIVATE kls_suitesparse)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(klu_compare PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(klu_compare PRIVATE
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

if(KLS_BUILD_SUBTREELU_COMPARE)
  if(NOT EXISTS "${SUBTREELU_ROOT}/include/subtree_lu.h")
    message(FATAL_ERROR
      "KLS_BUILD_SUBTREELU_COMPARE requires -DSUBTREELU_ROOT=/path/to/SubtreeLU")
  endif()
  enable_language(CXX)
  find_library(SUBTREELU_LIBRARY
    NAMES subtree_lu
    PATHS "${SUBTREELU_ROOT}/lib"
    NO_DEFAULT_PATH)
  if(NOT SUBTREELU_LIBRARY)
    message(FATAL_ERROR "Could not find libsubtree_lu under SUBTREELU_ROOT")
  endif()
  add_executable(subtreelu_compare
    bench/subtreelu_compare.cpp)
  target_include_directories(subtreelu_compare PRIVATE
    "${SUBTREELU_ROOT}/include")
  target_link_libraries(subtreelu_compare PRIVATE "${SUBTREELU_LIBRARY}" m)
  target_compile_features(subtreelu_compare PRIVATE cxx_std_17)
  target_compile_options(subtreelu_compare PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
endif()

if(KLS_BUILD_TESTS)
  enable_testing()
  if(Python3_Interpreter_FOUND)
    add_test(NAME kls_componentwise_differential
      COMMAND "${CMAKE_COMMAND}" -E env "CC=${CMAKE_C_COMPILER}"
              "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_componentwise_certificate.py")
    set_tests_properties(kls_componentwise_differential PROPERTIES
      RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
      add_test(
        NAME kls_accuracy_build_flags
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_accuracy_build_flags.py"
                "${CMAKE_C_COMPILER}")
      set_tests_properties(kls_accuracy_build_flags PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
    endif()
    add_test(
      NAME kls_slim_build_provenance
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_slim_build_provenance.py")
    add_test(
      NAME kls_license_boundary
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/scripts/audit_license_boundary.py"
              --root "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(
      NAME kls_generalization_manifests
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/scripts/audit_generalization_manifests.py"
              --development
              "${CMAKE_CURRENT_SOURCE_DIR}/bench/suitesparse_generalization_dev_manifest.txt"
              --holdout
              "${CMAKE_CURRENT_SOURCE_DIR}/bench/suitesparse_generalization_holdout_manifest.txt")
    add_test(
      NAME kls_large_generalization_manifests
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/scripts/audit_generalization_manifests.py"
              --development
              "${CMAKE_CURRENT_SOURCE_DIR}/bench/suitesparse_generalization_large_dev_manifest.txt"
              --holdout
              "${CMAKE_CURRENT_SOURCE_DIR}/bench/suitesparse_generalization_large_holdout_manifest.txt"
              --expected-count 12)
    if(KLS_BUILD_BENCHMARKS)
      add_test(NAME kls_bench_lifecycle_timing
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_bench_lifecycle_timing.py"
                --bench "$<TARGET_FILE:kls_bench>")
      set_tests_properties(kls_bench_lifecycle_timing PROPERTIES TIMEOUT 120)
    endif()
    if(KLS_BUILD_AUTOTUNER AND KLS_BUILD_BENCHMARKS)
      add_test(
        NAME kls_autotune_profile
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_autotune_profile.py"
                --autotune "$<TARGET_FILE:kls-autotune>"
                --bench "$<TARGET_FILE:kls_bench>"
                --source "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
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
  add_executable(kls_cancellation tests/kls_cancellation.c)
  target_link_libraries(kls_cancellation PRIVATE kls)
  add_test(NAME kls_cancellation COMMAND kls_cancellation)
  add_executable(kls_snb_profile tests/kls_snb_profile.c)
  target_link_libraries(kls_snb_profile PRIVATE kls)
  add_test(NAME kls_snb_profile COMMAND kls_snb_profile)
  set(KLS_PRIVATE_TESTS kls_compact_transpose kls_snb_product kls_snb_build kls_snb_calibration_layouts kls_egraph_batch_plan)
  set(KLS_ACCURACY_TESTS
    kls_accuracy_owner kls_accuracy_best kls_accuracy_gate kls_accuracy_recovery kls_componentwise_caller kls_application_managed kls_certificate_pool
    kls_accuracy_epochs kls_accuracy_frames kls_accuracy_rhs
    kls_accuracy_factor_epochs kls_accuracy_fused_epochs)
  list(APPEND KLS_PRIVATE_TESTS ${KLS_ACCURACY_TESTS})
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    list(APPEND KLS_PRIVATE_TESTS kls_snb_masked_bounds kls_trsm_masked_contract)
  endif()
  foreach(KLS_PRIVATE_TEST IN LISTS KLS_PRIVATE_TESTS)
    add_executable(${KLS_PRIVATE_TEST} tests/${KLS_PRIVATE_TEST}.c)
    target_link_libraries(${KLS_PRIVATE_TEST} PRIVATE kls)
    target_include_directories(${KLS_PRIVATE_TEST} PRIVATE
      $<TARGET_PROPERTY:kls,INCLUDE_DIRECTORIES>)
    target_compile_definitions(${KLS_PRIVATE_TEST} PRIVATE
      $<TARGET_PROPERTY:kls,COMPILE_DEFINITIONS>)
    target_compile_options(${KLS_PRIVATE_TEST} PRIVATE
      $<TARGET_PROPERTY:kls,COMPILE_OPTIONS>)
    add_test(NAME ${KLS_PRIVATE_TEST} COMMAND ${KLS_PRIVATE_TEST})
  endforeach()

  set_tests_properties(${KLS_ACCURACY_TESTS} PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  set_tests_properties(kls_egraph_batch_plan PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120)
  add_executable(kls_accuracy_certificate tests/kls_accuracy_certificate.c)
  add_executable(kls_componentwise tests/kls_componentwise.c)
  target_link_libraries(kls_componentwise PRIVATE kls)
  add_test(NAME kls_componentwise COMMAND kls_componentwise)
  set_tests_properties(kls_componentwise PROPERTIES RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  add_executable(kls_componentwise_faults tests/kls_componentwise_faults.c)
  target_link_libraries(kls_componentwise_faults PRIVATE m)
  add_test(NAME kls_componentwise_faults COMMAND kls_componentwise_faults)
  set_tests_properties(kls_componentwise_faults PROPERTIES RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  add_executable(kls_original_operator_cache tests/kls_original_operator_cache.c)
  target_link_libraries(kls_original_operator_cache PRIVATE m)
  add_test(NAME kls_original_operator_cache COMMAND kls_original_operator_cache)
  set_tests_properties(kls_original_operator_cache PROPERTIES RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  add_executable(kls_componentwise_parallel tests/kls_componentwise_parallel.c)
  target_link_libraries(kls_componentwise_parallel PRIVATE m Threads::Threads)
  if(TARGET OpenMP::OpenMP_C)
    target_link_libraries(kls_componentwise_parallel PRIVATE OpenMP::OpenMP_C)
  endif()
  add_test(NAME kls_componentwise_parallel COMMAND kls_componentwise_parallel)
  set_tests_properties(kls_componentwise_parallel PROPERTIES RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  add_executable(kls_accuracy_schedule tests/kls_accuracy_schedule.c)
  target_link_libraries(kls_accuracy_schedule PRIVATE m)
  add_test(NAME kls_accuracy_schedule COMMAND kls_accuracy_schedule)
  set_tests_properties(kls_accuracy_schedule PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  add_executable(kls_accuracy_sweep tests/kls_accuracy_sweep.c)
  add_executable(kls_interval_bounds tests/kls_interval_bounds.c)
  target_link_libraries(kls_interval_bounds PRIVATE m)
  add_test(NAME kls_interval_bounds COMMAND kls_interval_bounds)
  set_tests_properties(kls_interval_bounds PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  target_link_libraries(kls_accuracy_sweep PRIVATE m)
  add_test(NAME kls_accuracy_sweep COMMAND kls_accuracy_sweep)
  set_tests_properties(kls_accuracy_sweep PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")
  target_link_libraries(kls_accuracy_certificate PRIVATE kls)
  add_test(NAME kls_accuracy_certificate COMMAND kls_accuracy_certificate)
  set_tests_properties(kls_accuracy_certificate PROPERTIES
    RUN_SERIAL TRUE TIMEOUT 120 LABELS "accuracy")

  add_executable(kls_singleton_metadata tests/kls_singleton_metadata.c)
  target_link_libraries(kls_singleton_metadata PRIVATE kls_suitesparse)
  target_compile_options(kls_singleton_metadata PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
  add_test(NAME kls_singleton_metadata COMMAND kls_singleton_metadata)
  add_executable(kls_transformed_residual tests/kls_transformed_residual.c)
  target_link_libraries(kls_transformed_residual PRIVATE KLS::kls)
  if(KLS_MATH_LIBRARY)
    target_link_libraries(kls_transformed_residual PRIVATE "${KLS_MATH_LIBRARY}")
  endif()
  target_compile_options(kls_transformed_residual PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
  add_test(NAME kls_transformed_residual COMMAND kls_transformed_residual)
endif()

install(TARGETS kls
  EXPORT KLSTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

if(KLS_BUILD_AUTOTUNER)
  install(TARGETS kls-autotune RUNTIME DESTINATION bin)
endif()

install(DIRECTORY include/ DESTINATION include)
