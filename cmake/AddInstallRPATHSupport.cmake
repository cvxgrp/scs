# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
AddInstallRPATHSupport
----------------------

Add support to RPATH during installation to the project and the targets

.. command:: add_install_rpath_support

  Add support to RPATH during installation to the project::

  .. code-block:: cmake

    add_install_rpath_support([BIN_DIRS dir [dir]]
                              [LIB_DIRS dir [dir]]
                              [INSTALL_NAME_DIR [dir]]
                              [DEPENDS condition [condition]]
                              [USE_LINK_PATH])

  Normally (depending on the platform) when you install a shared
  library you can either specify its absolute path as the install name,
  or leave just the library name itself. In the former case the library
  will be correctly linked during run time by all executables and other
  shared libraries, but it must not change its install location. This
  is often the case for libraries installed in the system default
  library directory (e.g. ``/usr/lib``).
  In the latter case, instead, the library can be moved anywhere in the
  file system but at run time the dynamic linker must be able to find
  it. This is often accomplished by setting environmental variables
  (i.e. ``LD_LIBRARY_PATH`` on Linux).
  This procedure is usually not desirable for two main reasons:

  - by setting the variable you are changing the default behaviour
    of the dynamic linker thus potentially breaking executables (not as
    destructive as ``LD_PRELOAD``)
  - the variable will be used only by applications spawned by the shell
    and not by other processes.

  RPATH aims in solving the issues introduced by the second
  installation method. Using run-path dependent libraries you can
  create a directory structure containing executables and dependent
  libraries that users can relocate without breaking it.
  A run-path dependent library is a dependent library whose complete
  install name is not known when the library is created.
  Instead, the library specifies that the dynamic loader must resolve
  the library’s install name when it loads the executable that depends
  on the library. The executable or the other shared library will
  hardcode in the binary itself the additional search directories
  to be passed to the dynamic linker. This works great in conjunction
  with relative paths.
  This command will enable support to RPATH to your project.
  It will enable the following things:

   - If the project builds shared libraries it will generate a run-path
     enabled shared library, i.e. its install name will be resolved
     only at run time.
   - In all cases (building executables and/or shared libraries)
     dependent shared libraries with RPATH support will have their name
     resolved only at run time, by embedding the search path directly
     into the built binary.

  The command has the following parameters:

  Options:
   - ``USE_LINK_PATH``: if passed the command will automatically adds to
     the RPATH the path to all the dependent libraries.

  Arguments:
   - ``BIN_DIRS`` list of directories when the targets (executable and
     plugins) will be installed.
   - ``LIB_DIRS`` list of directories to be added to the RPATH. These
     directories will be added "relative" w.r.t. the ``BIN_DIRS`` and
     ``LIB_DIRS``.
   - ``INSTALL_NAME_DIR`` directory where the libraries will be installed.
     This variable will be used only if ``CMAKE_SKIP_RPATH`` or
     ``CMAKE_SKIP_INSTALL_RPATH`` is set to ``TRUE`` as it will set the
     ``INSTALL_NAME_DIR`` on all targets
   - ``DEPENDS`` list of conditions that should be ``TRUE`` to enable
     RPATH, for example ``FOO; NOT BAR``.

  Note: see https://gitlab.kitware.com/cmake/cmake/issues/16589 for further
  details.

.. command:: target_append_install_rpath

  Add extra paths to RPATH for a specific target::

  .. code-block:: cmake

    target_append_install_rpath(<target>
                                <INSTALL_DESTINATION destination>
                                [LIB_DIRS dir [dir]]
                                [DEPENDS condition [condition]])

  Arguments:
   - ``INSTALL_DESTINATION`` path where the target will be installed.
   - ``LIB_DIRS`` list of directories to be added to the RPATH. These
     directories will be added "relative" w.r.t. the ``INSTALL_DESTINATION``.
   - ``DEPENDS`` list of conditions that should be ``TRUE`` to enable
     RPATH, for example ``FOO; NOT BAR``.

#]=======================================================================]

include(CMakeParseArguments)


macro(__AddInstallRPATHSupport_GET_SYSTEM_LIB_DIRS _var)
  # Find system implicit lib directories
  set(${_var} ${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES})
  if(EXISTS "/etc/debian_version") # is this a debian system ?
    if(CMAKE_LIBRARY_ARCHITECTURE)
      list(APPEND ${_var} "/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
                          "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}")
    endif()
  endif()
endmacro()


macro(__AddInstallRPATHSupport_APPEND_RELATIVE_RPATH _var _bin_dir _lib_dir)
  file(RELATIVE_PATH _rel_path ${_bin_dir} ${_lib_dir})
  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    list(APPEND ${_var} "@loader_path/${_rel_path}")
  else()
    list(APPEND ${_var} "\$ORIGIN/${_rel_path}")
  endif()
endmacro()



function(ADD_INSTALL_RPATH_SUPPORT)

  set(_options USE_LINK_PATH)
  set(_oneValueArgs INSTALL_NAME_DIR)
  set(_multiValueArgs BIN_DIRS
                      LIB_DIRS
                      DEPENDS)

  cmake_parse_arguments(_ARS "${_options}"
                             "${_oneValueArgs}"
                             "${_multiValueArgs}"
                             "${ARGN}")

  # if either RPATH or INSTALL_RPATH is disabled
  # and the INSTALL_NAME_DIR variable is set, then hardcode the install name
  if(CMAKE_SKIP_RPATH OR CMAKE_SKIP_INSTALL_RPATH)
    if(DEFINED _ARS_INSTALL_NAME_DIR)
      set(CMAKE_INSTALL_NAME_DIR ${_ARS_INSTALL_NAME_DIR} PARENT_SCOPE)
    endif()
  endif()

  if (CMAKE_SKIP_RPATH OR (CMAKE_SKIP_INSTALL_RPATH AND CMAKE_SKIP_BUILD_RPATH))
    return()
  endif()


  set(_rpath_available 1)
  if(DEFINED _ARS_DEPENDS)
    foreach(_dep ${_ARS_DEPENDS})
      string(REGEX REPLACE " +" ";" _dep "${_dep}")
      if(NOT (${_dep}))
        set(_rpath_available 0)
      endif()
    endforeach()
  endif()

  if(_rpath_available)

    # Enable RPATH on OSX.
    set(CMAKE_MACOSX_RPATH TRUE PARENT_SCOPE)

    __AddInstallRPATHSupport_get_system_lib_dirs(_system_lib_dirs)

    # This is relative RPATH for libraries built in the same project
    foreach(lib_dir ${_ARS_LIB_DIRS})
      list(FIND _system_lib_dirs "${lib_dir}" isSystemDir)
      if("${isSystemDir}" STREQUAL "-1")
        foreach(bin_dir ${_ARS_LIB_DIRS} ${_ARS_BIN_DIRS})
          __AddInstallRPATHSupport_append_relative_rpath(CMAKE_INSTALL_RPATH ${bin_dir} ${lib_dir})
        endforeach()
      endif()
    endforeach()
    if(NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "")
      list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
    endif()
    set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} PARENT_SCOPE)

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ${_ARS_USE_LINK_PATH} PARENT_SCOPE)

  endif()

endfunction()


function(TARGET_APPEND_INSTALL_RPATH _target)
  set(_options )
  set(_oneValueArgs INSTALL_DESTINATION)
  set(_multiValueArgs LIB_DIRS
                      DEPENDS)

  if (CMAKE_SKIP_RPATH OR (CMAKE_SKIP_INSTALL_RPATH AND CMAKE_SKIP_BUILD_RPATH))
    return()
  endif()

  cmake_parse_arguments(_ARS "${_options}"
                             "${_oneValueArgs}"
                             "${_multiValueArgs}"
                             "${ARGN}")

  set(_rpath_available 1)
  if(DEFINED _ARS_DEPENDS)
    foreach(_dep ${_ARS_DEPENDS})
      string(REGEX REPLACE " +" ";" _dep "${_dep}")
      if(NOT (${_dep}))
        set(_rpath_available 0)
      endif()
    endforeach()
  endif()

  if(_rpath_available)

    __AddInstallRPATHSupport_get_system_lib_dirs(_system_lib_dirs)

    get_target_property(_current_rpath ${_target} INSTALL_RPATH)
    foreach(lib_dir ${_ARS_LIB_DIRS})
      list(FIND _system_lib_dirs "${lib_dir}" isSystemDir)
      if("${isSystemDir}" STREQUAL "-1")
        __AddInstallRPATHSupport_append_relative_rpath(_current_rpath ${_ARS_INSTALL_DESTINATION} ${lib_dir})
      endif()
    endforeach()
    set_target_properties(${_target} PROPERTIES INSTALL_RPATH "${_current_rpath}")
  endif()

endfunction()
