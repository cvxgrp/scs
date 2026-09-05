# Consumes an INSTALLED scs package (see .github/workflows/cmake_consume.yml).
# The installed package, not this fixture, must enforce its consumer floor.
cmake_minimum_required(VERSION 3.18)
project(scs_consume C)

set(SCS_COMPONENT scsdir CACHE STRING "scs library to consume (scsdir, scsmkl, ...)")

# CI: a consumer that resolved LAPACK at the wrong integer width before scs
# must be refused by a static package.
option(PREFIND_LAPACK_WRONG_WIDTH "find LAPACK with 4-byte integers first" OFF)
if(PREFIND_LAPACK_WRONG_WIDTH)
  set(BLA_SIZEOF_INTEGER 4)
  find_package(LAPACK REQUIRED)
endif()

find_package(scs CONFIG REQUIRED)
if(NOT scs::${SCS_COMPONENT} IN_LIST scs_LIBRARIES)
  message(FATAL_ERROR "scs_LIBRARIES=${scs_LIBRARIES} lacks scs::${SCS_COMPONENT}")
endif()
add_executable(consume consume.c)
target_link_libraries(consume PRIVATE scs::${SCS_COMPONENT})
