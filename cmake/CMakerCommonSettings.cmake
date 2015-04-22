set(CMakerCommonSettings_CALLED FALSE) # make sure this file is called once

if(NOT CMakerCommonSettings_CALLED)
  ## General settings ==========================================================
  # postfix, based on type
  set(CMAKE_DEBUG_POSTFIX "_d" CACHE STRING "postfix applied to debug build of libraries")
  set(CMAKE_RELEASE_POSTFIX "" CACHE STRING "postfix applied to release build of libraries")
  set(CMAKE_RELWITHDEBINFO_POSTFIX "_rd" CACHE STRING "postfix applied to release-with-debug-information libraries")
  set(CMAKE_MINSIZEREL_POSTFIX "_s" CACHE STRING "postfix applied to minimium-size-build libraries")

  # work out the postfix; required where we use OUTPUT_NAME
  if(CMAKE_BUILD_TYPE MATCHES Release)
    set(EXE_POSTFIX)
  elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    set(EXE_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(EXE_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
  elseif(CMAKE_BUILD_TYPE MATCHES MinSizeRel)
    set(EXE_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
  endif(CMAKE_BUILD_TYPE MATCHES Release)

  if(NOT THE_PROJECT_ROOT OR THE_PROJECT_ROOT STREQUAL "")
    cmaker_print_error("Please set THE_PROJECT_ROOT to your project root!!!")
  endif()

  if(NOT THE_LIB_RUNTIME_OUTPUT_DIRECTORY)
    set(THE_LIB_RUNTIME_OUTPUT_DIRECTORY ${THE_PROJECT_ROOT}/bin CACHE PATH "Target for the binaries")
    set(THE_LIB_LIBRARY_OUTPUT_DIRECTORY ${THE_PROJECT_ROOT}/lib CACHE PATH "Target for the libraries")
  endif()
  link_directories(${THE_LIB_LIBRARY_OUTPUT_DIRECTORY})

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${THE_LIB_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${THE_LIB_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${THE_LIB_RUNTIME_OUTPUT_DIRECTORY})

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${THE_LIB_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${THE_LIB_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${THE_LIB_RUNTIME_OUTPUT_DIRECTORY})
endif(NOT CMakerCommonSettings_CALLED)