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
## tribol_register_simple_library(
##    NAME     <name> 
##    DIR      <dir>                -- required path to root
##    LIBS     <list of libraries>  -- will use NAME if not specified
##    INCLUDE  <include_path        -- required path to include directory
##    DEPENDS_ON <deps...>          -- optional dependencies
##
## Sets up and registers a 'simple' library as a blt-registered target.
##
## A 'simple' library can be defined by specifying any of the following:
##   - an INCLUDE path to compile against
##   - a list of LIBS to link against (it is not an error if some are missing)
##   - a list of dependency targets 
##
## After calling this macro, the following variables will be defined
##   ${name}_FOUND        - Indicates whether the library was found
##   ${name}_INCLUDE_DIRS - The library include directories, if applicable
##   ${name}_LIBRARIES    - The paths to the link libraries, if applicable
##
## The macro also registers a blt-registered target named ${name}
##------------------------------------------------------------------------------
macro(tribol_register_simple_library)

  set(options)
  set(singleValueArgs NAME DIR INCLUDE )
  set(multiValueArgs LIBS DEPS)

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if ( NOT DEFINED arg_NAME )
    message( FATAL_ERROR "tribol_register_simple_library needs a NAME argument for the target!" )
  endif()

  if ( NOT DEFINED arg_DIR )
    message( FATAL_ERROR "tribol_register_simple_library needs a DIR argument for the directory!" )
  endif()
  
  if ( NOT EXISTS ${arg_DIR} )
    message( FATAL_ERROR "invalid dir for tribol_register_simple_library -- ${arg_DIR}!" )
  endif()

  string(TOUPPER ${arg_NAME} _uppercase_name)

  set(_req_vars)
  
  # setup the includes
  if ( DEFINED arg_INCLUDE AND EXISTS ${arg_INCLUDE} )
      set( ${_uppercase_name}_INCLUDE_DIRS ${arg_INCLUDE})
      blt_list_append(TO _req_vars ELEMENTS ${_uppercase_name}_INCLUDE_DIRS)
  endif()

  # setup the libraries
  if ( DEFINED arg_LIBS )
    blt_find_libraries( FOUND_LIBS ${_uppercase_name}_LIBRARIES
                        NAMES      ${arg_LIBS}
                        REQUIRED   FALSE
                        PATHS      ${arg_DIR}/lib)
    blt_list_append(TO _req_vars ELEMENTS ${_uppercase_name}_LIBRARIES)
  endif()                        

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(${_uppercase_name}  DEFAULT_MSG ${_req_vars})
  
  # mark library as found                                
  set(${_uppercase_name}_FOUND ${${_uppercase_name}_FOUND} CACHE BOOL "" FORCE)                                  

  mark_as_advanced( ${_uppercase_name}_FOUND
                    ${_req_vars} )

  blt_register_library(
    NAME ${arg_NAME}
    DEPENDS_ON ${arg_DEPS}
    INCLUDES ${${_uppercase_name}_INCLUDE_DIRS}
    LIBRARIES ${${_uppercase_name}_LIBRARIES}
    TREAT_INCLUDES_AS_SYSTEM TRUE )

endmacro(tribol_register_simple_library)


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
