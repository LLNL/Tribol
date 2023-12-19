# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

##------------------------------------------------------------------------------
## tribol_assert_path_exists( path )
##
## Checks if the specified path to a file or directory exits. If the path
## does not exist, CMake will throw a fatal error message and the configure
## step will fail.
##
##------------------------------------------------------------------------------
macro(tribol_assert_path_exists path )

  if ( NOT EXISTS ${path} )
    message( FATAL_ERROR "[${path}] does not exist!" )
  endif()

endmacro(tribol_assert_path_exists)


##------------------------------------------------------------------------------
## tribol_install( SOURCE_HEADERS    [ header1 [header2...] ]
##                 GENERATED_HEADERS [ h1 [h2...] ] )
##
## This macro installs the Tribol library, headers and CMake export targets.
##------------------------------------------------------------------------------
macro(tribol_install)

  set(options)
  set(singleValueArgs)
  set(multiValueArgs SOURCE_HEADERS GENERATED_HEADERS )

  ## parse arguments
  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if ( NOT arg_SOURCE_HEADERS )
    message( FATAL_ERROR "tribol_install() called without specifying SOURCE_HEADERS" )
  endif()

  install( TARGETS tribol
           EXPORT   tribol-targets
           LIBRARY  DESTINATION lib
           ARCHIVE  DESTINATION lib
           RUNTIME  DESTINATION bin
           INCLUDES DESTINATION include )

  ## mirror the directory structure on the install
  foreach( tribol_header_file ${arg_SOURCE_HEADERS} )
    get_filename_component( tribol_base_dir ${tribol_header_file} DIRECTORY)
    install( FILES ${tribol_header_file}
             DESTINATION "include/tribol/${tribol_base_dir}" )
  endforeach()

  ## copy over generated headers
  foreach( tribol_gen_header ${arg_GENERATED_HEADERS} )
     install( FILES ${tribol_gen_header}
              DESTINATION "include/tribol" )
  endforeach()

  install(EXPORT tribol-targets DESTINATION lib/cmake)

endmacro(tribol_install)

##------------------------------------------------------------------------------
## convert_to_native_escaped_file_path( path output )
##
## This macro converts a cmake path to a platform specific string literal
## usable in C++.  (For example, on windows C:/Path will be come C:\\Path)
##------------------------------------------------------------------------------

macro(convert_to_native_escaped_file_path path output)
    file(TO_NATIVE_PATH ${path} ${output})
    string(REPLACE "\\" "\\\\"  ${output} "${${output}}")
endmacro()

##------------------------------------------------------------------------------
## tribol_configure_file
##
## This macro is a thin wrapper over the builtin configure_file command.
## It has the same arguments/options as configure_file but introduces an
## intermediate file that is only copied to the target file if the target differs
## from the intermediate.
##------------------------------------------------------------------------------
macro(tribol_configure_file _source _target)
    set(_tmp_target ${_target}.tmp)
    configure_file(${_source} ${_tmp_target} ${ARGN})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_tmp_target} ${_target})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${_tmp_target})
endmacro(tribol_configure_file)
