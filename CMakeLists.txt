# CMakeLists.txt file for scs. This software may be modified and distributed
# under the terms of the MIT License

cmake_minimum_required(VERSION 3.5)

project(
  scs
  LANGUAGES C
  VERSION 3.2.6)

# Defines the CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_BINDIR and many other useful
# macros. See https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

# Control where libraries and executables are placed during the build. With the
# following settings executables are placed in <the top level of the build
# tree>/bin and libraries/archives in <top level of the build tree>/lib. This is
# particularly useful to run ctests on libraries built on Windows machines:
# tests, which are executables, are placed in the same folders of dlls, which
# are treated as executables as well, so that they can properly find the
# libraries to run. This is a because of missing RPATH on Windows.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# To build shared libraries in Windows, we set CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS
# to TRUE. See
# https://cmake.org/cmake/help/v3.4/variable/CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS.html
# See
# https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Under MSVC, we set CMAKE_DEBUG_POSTFIX to "d" to add a trailing "d" to library
# built in debug mode. In this Windows user can compile, build and install the
# library in both Release and Debug configuration avoiding naming clashes in the
# installation directories.
if(MSVC)
  set(CMAKE_DEBUG_POSTFIX "d")
endif()

# Build position independent code. Position Independent Code (PIC) is commonly
# used for shared libraries so that the same shared library code can be loaded
# in each program address space in a location where it will not overlap with any
# other uses of such memory. In particular, this option avoids problems
# occurring when a process wants to load more than one shared library at the
# same virtual address. Since shared libraries cannot predict where other shared
# libraries could be loaded, this is an unavoidable problem with the traditional
# shared library concept. Generating position-independent code is often the
# default behavior for most modern compilers. Moreover linking a static library
# that is not built with PIC from a shared library will fail on some
# compiler/architecture combinations. Further details on PIC can be found here:
# https://eli.thegreenplace.net/2011/11/03/position-independent-code-pic-in-shared-libraries/
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable C and C++ compiler extensions. C/CXX_EXTENSIONS are ON by default to
# allow the compilers to use extended variants of the C/CXX language. However,
# this could expose cross-platform bugs in user code or in the headers of
# third-party dependencies and thus it is strongly suggested to turn extensions
# off.
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_EXTENSIONS OFF)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options Shared/Dynamic or Static library?
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

# Encourage user to specify a build type (e.g. Release, Debug, etc.), otherwise
# set it to Release.
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
  endif()
endif()

# Build test related commands?
option(BUILD_TESTING "Create tests using CMake" OFF)
if(BUILD_TESTING)
  enable_testing()
endif()

# Add uninstall target After 'make install' can run 'make uninstall' to remove.
include(AddUninstallTarget)

# Some variables useful for sampling the building process. Note that the GPU
# profile is not compiled.
set(LINSYS linsys)
set(DIRSRC ${LINSYS}/cpu/direct)
set(INDIRSRC ${LINSYS}/cpu/indirect)
set(EXTERNAL ${LINSYS}/external)

# Options
# ----------------------------------------------
# Use floats instead of doubles
option(SFLOAT "Use single precision floats rather than doubles" OFF)
message(STATUS "Single precision floats (32bit) are ${SFLOAT}")

# Use long integers for indexing
option(DLONG "Use long integers (64bit) for indexing" OFF)
message(STATUS "Long integers (64bit) are ${DLONG}")

# Disable all printing
option(NO_PRINTING "Disables all printing" OFF)
message(STATUS "Printing is NOT ${NO_PRINTING}")

# Disable all read/write functionality
option(NO_READ_WRITE "Disables all read/write functionality" OFF)
message(STATUS "Read/write functionality is NOT ${NO_READ_WRITE}")

option(USE_LAPACK "Whether to use BLAS/LAPACK" ON)
message(STATUS "BLAS/LAPACK usage is ${USE_LAPACK}")

# Enable OpenMP support
option(USE_OPENMP "Compile with OpenMP support" OFF)
message(STATUS "OpenMP parallelization is ${USE_OPENMP}")

set(COMPILER_OPTS "-DCTRLC")
set(LAPACK_LINK_LIBRARIES "")

# Primitive types
if(SFLOAT)
  set(SCS_FLOAT_TYPE "float")
  set(COMPILER_OPTS "-DSFLOAT ${COMPILER_OPTS}")
else()
  set(SCS_FLOAT_TYPE "double")
endif()

if(DLONG)
  set(SCS_INT_TYPE "long long")
  set(COMPILER_OPTS "-DDLONG ${COMPILER_OPTS}")
else()
  set(SCS_INT_TYPE "int")
endif()

if(NO_PRINTING)
  set(COMPILER_OPTS "-DNO_PRINTING=1 ${COMPILER_OPTS}")
endif()

if(NO_READ_WRITE)
  set(COMPILER_OPTS "-DNO_READ_WRITE=1 ${COMPILER_OPTS}")
endif()

if (USE_LAPACK)
  set(COMPILER_OPTS "-DUSE_LAPACK ${COMPILER_OPTS}")
  list(APPEND LAPACK_LINK_LIBRARIES "blas" "lapack")
endif()

if(USE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif()

message(STATUS "COMPILER_OPTS = ${COMPILER_OPTS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# TODO this is a hack that overwrites the scs_types.h file, we should find a way
# to do this that doesn't pollute the master directory.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/scs_types.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/include/scs_types.h NEWLINE_STYLE LF)

# Public headers
set(${PROJECT_NAME}_PUBLIC_HDR include/scs_types.h include/scs.h)

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

# Common header files
set(${PROJECT_NAME}_HDR
    include/aa.h
    include/cones.h
    include/ctrlc.h
    include/glbopts.h
    include/linalg.h
    include/linsys.h
    include/normalize.h
    include/rw.h
    include/scs.h
    include/scs_blas.h
    include/scs_types.h
    include/scs_work.h
    include/util.h
    ${LINSYS}/csparse.h
    ${LINSYS}/scs_matrix.h)

# get all the c file in amd/external
file(GLOB ${PROJECT_NAME}_AMD_EXTERNAL_SRC ${EXTERNAL}/amd/*.c)

# get all the h file in amd/external
file(GLOB ${PROJECT_NAME}_AMD_EXTERNAL_HDR ${EXTERNAL}/amd/*.h)

# get all the c file in amd/external
file(GLOB ${PROJECT_NAME}_LDL_EXTERNAL_HDR ${EXTERNAL}/qdldl/*.h)

# ##############################################################################
# Direct method. Here we compile the direct method library.
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
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scs>")

# Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, ECPs, and PCPs
target_compile_definitions(${${PROJECT_NAME}_DIRECT} PRIVATE ${COMPILER_OPTS})

# The library depends on math (m) and (optionally) blas and lapack
target_link_libraries(${${PROJECT_NAME}_DIRECT} PRIVATE m ${LAPACK_LINK_LIBRARIES})

# Set some properties
set_target_properties(
  ${${PROJECT_NAME}_DIRECT}
  PROPERTIES VERSION ${scs_VERSION} PUBLIC_HEADER
                                    "${${PROJECT_NAME}_PUBLIC_HDR}")

add_library(scs::${${PROJECT_NAME}_DIRECT} ALIAS ${${PROJECT_NAME}_DIRECT})

# Install the library
install(
  TARGETS ${${PROJECT_NAME}_DIRECT}
  EXPORT ${PROJECT_NAME}
  COMPONENT runtime
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

# Add the direct method to the set of the compiled targets
set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_DIRECT})

# ##############################################################################

# Indirect method. Here we compile the indirect method library.
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
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scs>")

# Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, ECPs, and PCPs
target_compile_definitions(${${PROJECT_NAME}_INDIRECT} PRIVATE ${COMPILER_OPTS}
                                                               -DINDIRECT)

# The library depends on math (m) and (optionally) blas and lapack
target_link_libraries(${${PROJECT_NAME}_INDIRECT} PUBLIC m ${LAPACK_LINK_LIBRARIES})

# Set some properties
set_target_properties(
  ${${PROJECT_NAME}_INDIRECT}
  PROPERTIES VERSION ${scs_VERSION} PUBLIC_HEADER
                                    "${${PROJECT_NAME}_PUBLIC_HDR}")

add_library(scs::${${PROJECT_NAME}_INDIRECT} ALIAS ${${PROJECT_NAME}_INDIRECT})

# Install the library
install(
  TARGETS ${${PROJECT_NAME}_INDIRECT}
  EXPORT ${PROJECT_NAME}
  COMPONENT runtime
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

# Add the indirect method to the set of the compiled targets
set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_INDIRECT})

# ##############################################################################

# If MKLROOT is defined, then we install the MKL version.
if(DEFINED ENV{MKLROOT})

  message(STATUS "MKLROOT set to $ENV{MKLROOT}")
  message(STATUS "Will install SCS-MKL (libscsmkl).")

  set(MKLSRC ${LINSYS}/mkl/direct)

  # Here we compile the direct MKL pardiso library
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
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scs>")

  target_compile_definitions(${${PROJECT_NAME}_MKL} PRIVATE ${COMPILER_OPTS})
  # See:
  # https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
  # This is probably not correct for other systems. TODO: make SCS-MKL
  # work for all combinations of platform / compiler / threading options.
  target_link_options(${${PROJECT_NAME}_MKL} PRIVATE "LINKER:--no-as-needed")
  target_link_directories(${${PROJECT_NAME}_MKL} PRIVATE $ENV{MKLROOT}/lib)
  target_link_libraries(
    ${${PROJECT_NAME}_MKL}
    PRIVATE m
            mkl_rt
            mkl_gnu_thread
            mkl_core
            gomp
            pthread
            dl)

  # Set some properties
  set_target_properties(
    ${${PROJECT_NAME}_MKL}
    PROPERTIES VERSION ${scs_VERSION} PUBLIC_HEADER
                                      "${${PROJECT_NAME}_PUBLIC_HDR}")

  add_library(scs::${${PROJECT_NAME}_MKL} ALIAS ${${PROJECT_NAME}_MKL})

  # Install the library
  install(
    TARGETS ${${PROJECT_NAME}_MKL}
    EXPORT ${PROJECT_NAME}
    COMPONENT runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/scs")

  # Add the mkl method to the set of the compiled targets
  set_property(GLOBAL APPEND PROPERTY SCS_TARGETS ${${PROJECT_NAME}_MKL})
else()
  message(
    STATUS "MKLROOT environment variable is undefined, skipping SCS-MKL install"
  )
endif()

# ##############################################################################

# Install the .cmake file. It is possible to link scs o consumer libraries.
include(InstallBasicPackageFiles)
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
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Add the tests
if(BUILD_TESTING)
  add_executable(run_tests_direct test/run_tests.c)
  target_compile_definitions(run_tests_direct PRIVATE ${COMPILER_OPTS})
  target_link_libraries(run_tests_direct PRIVATE scs::scsdir)
  target_include_directories(
    run_tests_direct
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test;${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  add_test(
    NAME run_tests_direct
    COMMAND run_tests_direct
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(run_tests_indirect test/run_tests.c)
  target_compile_definitions(run_tests_indirect PRIVATE ${COMPILER_OPTS})
  target_link_libraries(run_tests_indirect PRIVATE scs::scsindir)
  target_include_directories(
    run_tests_indirect
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test;${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
  )

  add_test(
    NAME run_tests_indirect
    COMMAND run_tests_indirect
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  if(DEFINED ENV{MKLROOT})
    add_executable(run_tests_mkl test/run_tests.c)
    target_compile_definitions(run_tests_mkl PRIVATE ${COMPILER_OPTS})
    target_link_libraries(run_tests_mkl PRIVATE scs::scsmkl)
    target_include_directories(
      run_tests_mkl
      PRIVATE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test;${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_SOURCE_DIR}/${LINSYS}>"
    )

    add_test(
      NAME run_tests_mkl
      COMMAND run_tests_mkl
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

endif()
