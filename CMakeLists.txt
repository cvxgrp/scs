# CMakeLists.txt file for scs. This software may be modified and distributed
# under the terms of the MIT License

# 3.22: FindLAPACK honors BLA_SIZEOF_INTEGER (the LP64/ILP64 contract) and
# provides the LAPACK::LAPACK target that keeps static packages relocatable.
cmake_minimum_required(VERSION 3.22)

project(
  scs
  LANGUAGES C
  VERSION 3.3.1)

# ABI version, the loader identity of the shared libraries (libscsdir.so.3.3):
# bump on ABI breaks, never on patch releases. scs.mk derives the same value
# from SCS_VERSION in include/glbopts.h.
set(SCS_ABI_VERSION "${scs_VERSION_MAJOR}.${scs_VERSION_MINOR}")

# Standard GNU install directories (lib, bin, include, etc.)
include(GNUInstallDirs)

# Place executables in bin/ and libraries in lib/ within the build tree.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Export all symbols on Windows (avoids __declspec(dllexport) annotations).
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Append "d" suffix to debug libraries under MSVC to avoid name clashes.
if(MSVC)
  set(CMAKE_DEBUG_POSTFIX "d")
endif()

# Position-independent code is required for shared libraries.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable compiler extensions for strict standard compliance.
# Strict C99, the oldest standard a shipping consumer compiles scs with (the
# MATLAB binding uses -std=c99): GNU-only constructs fail here first.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
# Strict C99 hides POSIX declarations on glibc and illumos (sigaction,
# clock_gettime, struct timespec) unless a feature-test macro asks for them.
# A compile definition reaches every translation unit whatever its include
# order, which a header cannot guarantee.
if(CMAKE_SYSTEM_NAME MATCHES "Linux|SunOS")
  add_compile_definitions(_POSIX_C_SOURCE=200809L)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)

# Enable RPATH support for installed binaries and libraries
include(AddInstallRPATHSupport)
add_install_rpath_support(
  BIN_DIRS
  "${CMAKE_INSTALL_FULL_BINDIR}"
  LIB_DIRS
  "${CMAKE_INSTALL_FULL_LIBDIR}"
  INSTALL_NAME_DIR
  "${CMAKE_INSTALL_FULL_LIBDIR}"
  USE_LINK_PATH)

# Default to Release build type if none specified.
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
  endif()
endif()

option(BUILD_TESTING "Create tests using CMake" OFF)
if(BUILD_TESTING)
  enable_testing()
endif()

# Support 'make uninstall' after install.
include(AddUninstallTarget)

# Linear system solver directories.
set(LINSYS linsys)
set(DIRSRC ${LINSYS}/cpu/direct)
set(INDIRSRC ${LINSYS}/cpu/indirect)
set(EXTERNAL ${LINSYS}/external)

# Build options
option(SFLOAT "Use single precision floats rather than doubles" OFF)
message(STATUS "Single precision floats (32bit) are ${SFLOAT}")

option(DLONG "Use long integers (64bit) for indexing" OFF)
message(STATUS "Long integers (64bit) are ${DLONG}")

option(NO_PRINTING "Disables all printing" OFF)
message(STATUS "Printing is NOT ${NO_PRINTING}")

option(NO_READ_WRITE "Disables all read/write functionality" OFF)
message(STATUS "Read/write functionality is NOT ${NO_READ_WRITE}")

option(USE_LAPACK "Whether to use BLAS/LAPACK" ON)
message(STATUS "BLAS/LAPACK usage is ${USE_LAPACK}")

option(BLAS64 "Use 64-bit integers for BLAS/LAPACK" OFF)
message(STATUS "64-bit BLAS integers are ${BLAS64}")

option(USE_SPECTRAL_CONES "Whether to use spectral cones (requires LAPACK)" OFF)
if(USE_SPECTRAL_CONES AND NOT USE_LAPACK)
  message(FATAL_ERROR "USE_SPECTRAL_CONES requires USE_LAPACK to be ON")
endif()

option(USE_OPENMP "Compile with OpenMP support" OFF)
message(STATUS "OpenMP parallelization is ${USE_OPENMP}")

option(USE_CUDSS "Compile with cuDSS support" OFF)
message(STATUS "cuDSS support is ${USE_CUDSS}")

if(APPLE)
  option(USE_APPLE_ACCELERATE "Use Apple Accelerate framework for sparse LDLt" OFF)
  message(STATUS "Apple Accelerate support is ${USE_APPLE_ACCELERATE}")
endif()

set(COMPILER_OPTS CTRLC)
set(LAPACK_LINK_LIBRARIES "")
set(SCS_TEST_LINK_LIBRARIES "")
if(NOT WIN32)
  list(APPEND SCS_TEST_LINK_LIBRARIES m)
endif()

# Primitive types
if(SFLOAT)
  set(SCS_FLOAT_TYPE "float")
  list(APPEND COMPILER_OPTS SFLOAT)
else()
  set(SCS_FLOAT_TYPE "double")
endif()

if(DLONG)
  set(SCS_INT_TYPE "long long")
  list(APPEND COMPILER_OPTS DLONG)
else()
  set(SCS_INT_TYPE "int")
endif()

if(NO_PRINTING)
  list(APPEND COMPILER_OPTS NO_PRINTING=1)
endif()

if(NO_READ_WRITE)
  list(APPEND COMPILER_OPTS NO_READ_WRITE=1)
endif()

if(USE_LAPACK)
  list(APPEND COMPILER_OPTS USE_LAPACK)
  if(BLAS64)
    list(APPEND COMPILER_OPTS BLAS64)
    set(BLA_SIZEOF_INTEGER 8)
  else()
    set(BLA_SIZEOF_INTEGER 4)
  endif()
  if(DEFINED ENV{MKLROOT} AND NOT DEFINED BLA_VENDOR)
    if(BLAS64)
      set(BLA_VENDOR Intel10_64ilp)
    else()
      set(BLA_VENDOR Intel10_64lp)
    endif()
  endif()
  find_package(LAPACK)
  if(LAPACK_FOUND)
    # Link the imported target, not resolved paths: a static library exports
    # its private link line, and absolute provider paths would make the
    # installed package usable only on this machine.
    set(LAPACK_LINK_LIBRARIES LAPACK::LAPACK)
  else()
    if(BLAS64)
      message(FATAL_ERROR
              "BLAS64=ON requires a discoverable ILP64 BLAS/LAPACK library. "
              "CMake did not find one. Set BLA_VENDOR explicitly or configure "
              "with -DBLAS64=OFF.")
    endif()
    # Fallback to linking by name (works on most Linux/Mac systems).
    set(LAPACK_LINK_LIBRARIES blas lapack)
  endif()
endif()

find_package(Threads REQUIRED)

if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_C_FOUND)
    # Link target used below via target_link_libraries.
  else()
    message(WARNING "OpenMP requested but not found by CMake")
  endif()
endif()

if(USE_SPECTRAL_CONES)
  list(APPEND COMPILER_OPTS USE_SPECTRAL_CONES)
endif()

option(USE_SPECTRAL_TIMING_FLAG "Enable timing of spectral cone projections" OFF)
if(USE_SPECTRAL_TIMING_FLAG)
  list(APPEND COMPILER_OPTS SPECTRAL_TIMING_FLAG)
endif()

option(USE_SPECTRAL_DEBUG "Enable spectral cone debug output" OFF)
if(USE_SPECTRAL_DEBUG)
  list(APPEND COMPILER_OPTS SPECTRAL_DEBUG)
endif()

message(STATUS "COMPILER_OPTS = ${COMPILER_OPTS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Generate a build-tree include overlay so CMake configuration never modifies
# the checked-in Make-compatible headers.
set(SCS_BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
set(SCS_GENERATED_TYPES_HDR "${SCS_BUILD_INCLUDE_DIR}/scs_types.h")
file(MAKE_DIRECTORY "${SCS_BUILD_INCLUDE_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/scs_types.h.in
               ${SCS_GENERATED_TYPES_HDR} NEWLINE_STYLE LF)

set(SCS_CONFIGURED_INCLUDE_HDR
    include/aa.h
    include/aa_stats.h
    include/cones.h
    include/ctrlc.h
    include/glbopts.h
    include/linalg.h
    include/linsys.h
    include/normalize.h
    include/rw.h
    include/scs.h
    include/scs_blas.h
    include/scs_work.h
    include/util.h
    include/util_spectral_cones.h)

set(SCS_BUILD_INCLUDE_HDR "")
foreach(_scs_header ${SCS_CONFIGURED_INCLUDE_HDR})
  get_filename_component(_scs_header_name "${_scs_header}" NAME)
  set(_scs_build_header "${SCS_BUILD_INCLUDE_DIR}/${_scs_header_name}")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${_scs_header}"
                 "${_scs_build_header}" COPYONLY)
  list(APPEND SCS_BUILD_INCLUDE_HDR "${_scs_build_header}")
endforeach()

# Public headers
set(${PROJECT_NAME}_PUBLIC_HDR
    "${SCS_BUILD_INCLUDE_DIR}/aa_stats.h"
    "${SCS_GENERATED_TYPES_HDR}"
    "${SCS_BUILD_INCLUDE_DIR}/scs.h")
set(SCS_PUBLIC_INCLUDE_DIRS
    "$<BUILD_INTERFACE:${SCS_BUILD_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scs>")

# Common source files
set(${PROJECT_NAME}_SRC
    src/aa.c
    src/cones.c
    src/ctrlc.c
    src/exp_cone.c
    src/linalg.c
    src/normalize.c
    src/rw.c
    src/scs.c
    src/scs_version.c
    src/util.c
    ${LINSYS}/csparse.c
    ${LINSYS}/scs_matrix.c)

# Spectral source files
set(SPECTRAL_SRC
    src/spectral_cones/logdeterminant/log_cone_Newton.c
    src/spectral_cones/logdeterminant/log_cone_IPM.c
    src/spectral_cones/logdeterminant/log_cone_wrapper.c
    src/spectral_cones/logdeterminant/logdet_cone.c
    src/spectral_cones/nuclear/ell1_cone.c
    src/spectral_cones/nuclear/nuclear_cone.c
    src/spectral_cones/sum-largest/sum_largest_cone.c
    src/spectral_cones/sum-largest/sum_largest_eval_cone.c
    src/spectral_cones/util_spectral_cones.c)

# Conditionally add spectral files
if(USE_SPECTRAL_CONES)
  list(APPEND ${PROJECT_NAME}_SRC ${SPECTRAL_SRC})
endif()

# Common header files
set(${PROJECT_NAME}_HDR
    ${SCS_BUILD_INCLUDE_HDR}
    ${SCS_GENERATED_TYPES_HDR}
    ${LINSYS}/csparse.h
    ${LINSYS}/scs_matrix.h)

# Spectral header files
set(SPECTRAL_HDR "${SCS_BUILD_INCLUDE_DIR}/util_spectral_cones.h")

if(USE_SPECTRAL_CONES)
  list(APPEND ${PROJECT_NAME}_HDR ${SPECTRAL_HDR})
endif()

# External vendored dependencies (AMD ordering, QDLDL factorization).
file(GLOB ${PROJECT_NAME}_AMD_EXTERNAL_SRC ${EXTERNAL}/amd/*.c)
file(GLOB ${PROJECT_NAME}_AMD_EXTERNAL_HDR ${EXTERNAL}/amd/*.h)
file(GLOB ${PROJECT_NAME}_LDL_EXTERNAL_HDR ${EXTERNAL}/qdldl/*.h)

# Direct method library.
set(${PROJECT_NAME}_DIRECT ${PROJECT_NAME}dir)
add_library(
  ${${PROJECT_NAME}_DIRECT}
  ${${PROJECT_NAME}_HDR}
  ${${PROJECT_NAME}_SRC}
  ${DIRSRC}/private.c
  ${DIRSRC}/private.h
  ${EXTERNAL}/qdldl/qdldl.c
  ${${PROJECT_NAME}_AMD_EXTERNAL_SRC}
  ${${PROJECT_NAME}_AMD_EXTERNAL_HDR}
  ${${PROJECT_NAME}_LDL_EXTERNAL_HDR})

target_include_directories(
  ${${PROJECT_NAME}_DIRECT}
  PRIVATE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${EXTERNAL}/amd;${CMAKE_CURRENT_SOURCE_DIR}/${EXTERNAL}/qdldl;${CMAKE_CURRENT_SOURCE_DIR}/${DIRSRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
)

target_include_directories(
  ${${PROJECT_NAME}_DIRECT}
  PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

target_compile_definitions(${${PROJECT_NAME}_DIRECT} PRIVATE ${COMPILER_OPTS})

target_link_libraries(${${PROJECT_NAME}_DIRECT}
                      PRIVATE ${LAPACK_LINK_LIBRARIES}
                              $<$<NOT:$<PLATFORM_ID:Windows>>:m>
                              $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>
                              Threads::Threads)

set_target_properties(
  ${${PROJECT_NAME}_DIRECT}
  PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                    "${${PROJECT_NAME}_PUBLIC_HDR}")

add_library(scs::${${PROJECT_NAME}_DIRECT} ALIAS ${${PROJECT_NAME}_DIRECT})

install(
  TARGETS ${${PROJECT_NAME}_DIRECT}
  EXPORT ${PROJECT_NAME}
  COMPONENT runtime
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_DIRECT})

# Indirect method library.
set(${PROJECT_NAME}_INDIRECT ${PROJECT_NAME}indir)
add_library(
  ${${PROJECT_NAME}_INDIRECT}
  ${${PROJECT_NAME}_HDR} ${${PROJECT_NAME}_SRC} ${INDIRSRC}/private.c
  ${INDIRSRC}/private.h ${${${PROJECT_NAME}_INDIRECT}_HDR})

target_include_directories(
  ${${PROJECT_NAME}_INDIRECT}
  PRIVATE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INDIRSRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
)

target_include_directories(
  ${${PROJECT_NAME}_INDIRECT}
  PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

target_compile_definitions(${${PROJECT_NAME}_INDIRECT} PRIVATE ${COMPILER_OPTS}
                                                               -DINDIRECT)

target_link_libraries(${${PROJECT_NAME}_INDIRECT}
                      PRIVATE $<$<NOT:$<PLATFORM_ID:Windows>>:m>
                              ${LAPACK_LINK_LIBRARIES}
                              $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>
                              Threads::Threads)

set_target_properties(
  ${${PROJECT_NAME}_INDIRECT}
  PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                    "${${PROJECT_NAME}_PUBLIC_HDR}")

add_library(scs::${${PROJECT_NAME}_INDIRECT} ALIAS ${${PROJECT_NAME}_INDIRECT})

install(
  TARGETS ${${PROJECT_NAME}_INDIRECT}
  EXPORT ${PROJECT_NAME}
  COMPONENT runtime
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_INDIRECT})

# Dense direct method library (optional, requires LAPACK).
if(USE_LAPACK)
  set(DENSESRC ${LINSYS}/cpu/dense)

  set(${PROJECT_NAME}_DENSE ${PROJECT_NAME}dense)
  add_library(
    ${${PROJECT_NAME}_DENSE}
    ${${PROJECT_NAME}_HDR} ${${PROJECT_NAME}_SRC} ${DENSESRC}/private.c
    ${DENSESRC}/private.h)

  target_include_directories(
    ${${PROJECT_NAME}_DENSE}
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${DENSESRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  target_include_directories(
    ${${PROJECT_NAME}_DENSE}
    PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

  target_compile_definitions(${${PROJECT_NAME}_DENSE} PRIVATE ${COMPILER_OPTS})

  target_link_libraries(${${PROJECT_NAME}_DENSE}
                        PRIVATE ${LAPACK_LINK_LIBRARIES}
                                $<$<NOT:$<PLATFORM_ID:Windows>>:m>
                                $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>
                                Threads::Threads)

  set_target_properties(
    ${${PROJECT_NAME}_DENSE}
    PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                      "${${PROJECT_NAME}_PUBLIC_HDR}")

  add_library(scs::${${PROJECT_NAME}_DENSE} ALIAS ${${PROJECT_NAME}_DENSE})

  install(
    TARGETS ${${PROJECT_NAME}_DENSE}
    EXPORT ${PROJECT_NAME}
    COMPONENT runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

  set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_DENSE})
endif()

# MKL Pardiso library (optional, requires MKLROOT environment variable).
if(DEFINED ENV{MKLROOT})

  message(STATUS "MKLROOT set to $ENV{MKLROOT}")
  message(STATUS "Will install SCS-MKL (libscsmkl).")

  if(BLAS64 AND NOT DLONG)
    message(FATAL_ERROR
            "MKL PARDISO with BLAS64 requires DLONG=ON so that the "
            "pardiso_64 entry point matches the ILP64 BLAS interface.")
  endif()

  set(MKLSRC ${LINSYS}/mkl/direct)

  # MKL Pardiso direct solver library
  set(${PROJECT_NAME}_MKL ${PROJECT_NAME}mkl)
  add_library(
    ${${PROJECT_NAME}_MKL}
    ${${PROJECT_NAME}_HDR} ${${PROJECT_NAME}_SRC} ${MKLSRC}/private.c
    ${MKLSRC}/private.h ${${${PROJECT_NAME}_MKL}_HDR})

  target_include_directories(
    ${${PROJECT_NAME}_MKL}
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${MKLSRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  target_include_directories(
    ${${PROJECT_NAME}_MKL}
    PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

  target_compile_definitions(
    ${${PROJECT_NAME}_MKL}
    PRIVATE ${COMPILER_OPTS} SCS_MKL)
  # Link against mkl_rt (the single dynamic library) which dispatches to the
  # correct interface layer at runtime. The BLAS integer width (LP64 vs ILP64)
  # is set at startup by the MKL-specific path in src/scs.c, controlled by the
  # BLAS64 compile flag.
  #
  # Note: The PARDISO integer width (pardiso vs pardiso_64) is controlled by
  # DLONG, not by the MKL interface layer. See linsys/mkl/direct/private.c.
  # oneMKL's own package config (MKLConfig.cmake, oneMKL >= 2021.3) provides
  # the MKL::MKL target, which keeps a static export relocatable; SDL mode
  # links the single dynamic library (mkl_rt) whose MKL_Set_Interface_Layer
  # the runtime check uses. --no-as-needed keeps mkl_rt as DT_NEEDED for that
  # check; it is a GNU/ELF option, hence Linux only.
  set(MKL_LINK sdl)
  find_package(MKL CONFIG REQUIRED PATHS "$ENV{MKLROOT}")
  set(SCS_MKL_USES_CONFIG TRUE)
  target_link_options(${${PROJECT_NAME}_MKL} PRIVATE
                      "$<$<PLATFORM_ID:Linux>:LINKER:--no-as-needed>")
  target_link_libraries(
    ${${PROJECT_NAME}_MKL}
    PRIVATE $<$<NOT:$<PLATFORM_ID:Windows>>:m> MKL::MKL Threads::Threads
            ${CMAKE_DL_LIBS})

  set_target_properties(
    ${${PROJECT_NAME}_MKL}
    PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                      "${${PROJECT_NAME}_PUBLIC_HDR}")

  add_library(scs::${${PROJECT_NAME}_MKL} ALIAS ${${PROJECT_NAME}_MKL})

  install(
    TARGETS ${${PROJECT_NAME}_MKL}
    EXPORT ${PROJECT_NAME}
    COMPONENT runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

  set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_MKL})
else()
  message(
    STATUS "MKLROOT environment variable is undefined, skipping SCS-MKL install"
  )
endif()

# cuDSS GPU solver library (optional, requires CUDA toolkit and cuDSS).
if(USE_CUDSS)
  find_package(CUDAToolkit REQUIRED)
  find_package(cudss REQUIRED)

  # Force 32-bit integers for cuDSS compatibility
  if(DLONG)
    message(
      FATAL_ERROR
        "cuDSS requires 32-bit integers. Set DLONG=OFF or do not use -DDLONG=ON"
    )
  endif()

  message(STATUS "Will install SCS-cuDSS (libscscudss).")

  set(CUDSSSRC ${LINSYS}/cudss/direct)

  set(${PROJECT_NAME}_CUDSS ${PROJECT_NAME}cudss)
  add_library(
    ${${PROJECT_NAME}_CUDSS} ${${PROJECT_NAME}_HDR} ${${PROJECT_NAME}_SRC}
                             ${CUDSSSRC}/private.c ${CUDSSSRC}/private.h)

  target_include_directories(
    ${${PROJECT_NAME}_CUDSS}
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${CUDSSSRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  target_include_directories(
    ${${PROJECT_NAME}_CUDSS}
    PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

  target_compile_definitions(${${PROJECT_NAME}_CUDSS} PRIVATE ${COMPILER_OPTS})

  target_link_libraries(
    ${${PROJECT_NAME}_CUDSS}
    PRIVATE CUDA::cudart cudss $<$<NOT:$<PLATFORM_ID:Windows>>:m>
            ${LAPACK_LINK_LIBRARIES}
            Threads::Threads)

  set_target_properties(
    ${${PROJECT_NAME}_CUDSS}
    PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                      "${${PROJECT_NAME}_PUBLIC_HDR}")

  add_library(scs::${${PROJECT_NAME}_CUDSS} ALIAS ${${PROJECT_NAME}_CUDSS})

  install(
    TARGETS ${${PROJECT_NAME}_CUDSS}
    EXPORT ${PROJECT_NAME}
    COMPONENT runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

  set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_CUDSS})
endif()

# Apple Accelerate library (optional, macOS only).
if(APPLE AND USE_APPLE_ACCELERATE)
  if(DLONG)
    message(
      FATAL_ERROR
        "Apple Accelerate requires 32-bit integers. Set DLONG=OFF or do not use -DDLONG=ON"
    )
  endif()

  message(STATUS "Will install SCS-Accelerate (libscsaccel).")

  set(ACCELSRC ${LINSYS}/accelerate/direct)

  set(${PROJECT_NAME}_ACCEL ${PROJECT_NAME}accel)
  add_library(
    ${${PROJECT_NAME}_ACCEL}
    ${${PROJECT_NAME}_HDR} ${${PROJECT_NAME}_SRC} ${ACCELSRC}/private.c
    ${ACCELSRC}/private.h)

  target_include_directories(
    ${${PROJECT_NAME}_ACCEL}
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${ACCELSRC};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  target_include_directories(
    ${${PROJECT_NAME}_ACCEL}
    PUBLIC ${SCS_PUBLIC_INCLUDE_DIRS})

  target_compile_definitions(${${PROJECT_NAME}_ACCEL} PRIVATE ${COMPILER_OPTS})

  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(
    ${${PROJECT_NAME}_ACCEL}
    PRIVATE ${ACCELERATE_FRAMEWORK}
            $<$<BOOL:${OpenMP_C_FOUND}>:OpenMP::OpenMP_C>
            Threads::Threads)

  set_target_properties(
    ${${PROJECT_NAME}_ACCEL}
    PROPERTIES VERSION ${scs_VERSION} SOVERSION ${SCS_ABI_VERSION}
             PUBLIC_HEADER
                                      "${${PROJECT_NAME}_PUBLIC_HDR}")

  add_library(scs::${${PROJECT_NAME}_ACCEL} ALIAS ${${PROJECT_NAME}_ACCEL})

  install(
    TARGETS ${${PROJECT_NAME}_ACCEL}
    EXPORT ${PROJECT_NAME}
    COMPONENT runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

  set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_ACCEL})
endif()

# Install CMake package config files for find_package(scs) support.
include(InstallBasicPackageFiles)
# The config template re-finds a static library's private dependencies at the
# LAPACK integer width scs was built with; see cmake/scsConfig.cmake.in.
set(SCS_LAPACK_INT_BYTES "")
if(USE_LAPACK AND LAPACK_FOUND)
  set(SCS_LAPACK_INT_BYTES "${BLA_SIZEOF_INTEGER}")
endif()
install_basic_package_files(
  ${PROJECT_NAME}
  NAMESPACE
  scs::
  VERSION
  ${${PROJECT_NAME}_VERSION}
  TARGETS_PROPERTY
  SCS_TARGETS
  COMPATIBILITY
  SameMajorVersion
  VARS_PREFIX
  ${PROJECT_NAME}
  CONFIG_TEMPLATE
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/scsConfig.cmake.in"
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Tests
if(BUILD_TESTING)
  # Helper to avoid repeating test boilerplate for each backend.
  function(scs_add_test TEST_NAME LINK_TARGET)
    add_executable(${TEST_NAME} test/run_tests.c)
    target_compile_definitions(${TEST_NAME} PRIVATE ${COMPILER_OPTS})
    target_link_libraries(${TEST_NAME} PRIVATE ${LINK_TARGET}
                                               ${SCS_TEST_LINK_LIBRARIES})
    target_include_directories(
      ${TEST_NAME}
      PRIVATE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test;${SCS_BUILD_INCLUDE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
    )
    add_test(
      NAME ${TEST_NAME}
      COMMAND ${TEST_NAME}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endfunction()

  scs_add_test(run_tests_direct scs::scsdir)
  scs_add_test(run_tests_indirect scs::scsindir)

  if(USE_LAPACK)
    scs_add_test(run_tests_dense scs::scsdense)
  endif()

  if(DEFINED ENV{MKLROOT})
    scs_add_test(run_tests_mkl scs::scsmkl)
  endif()

  if(USE_CUDSS)
    scs_add_test(run_tests_cudss scs::scscudss)
  endif()

  if(APPLE AND USE_APPLE_ACCELERATE)
    scs_add_test(run_tests_accelerate scs::scsaccel)
  endif()
endif()
