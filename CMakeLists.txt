cmake_minimum_required(VERSION 2.8.9)

set(Closure_PROJECT_NAME closure_and_apps)
project(${Closure_PROJECT_NAME})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)  # for organizing code to folders

# Look for our CMake modules paths
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
                              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules
                              ${CMAKE_CURRENT_BINARY_DIR})
include(ClosureUtil)  # include useful macros
closure_print_status("Closure cmake modules: ${Closure_COLOR_BLUE}${CMAKE_MODULE_PATH}${Closure_COLOR_RESET}")

set(Closure_WORLD_VERSION 1)
set(Closure_MAJOR_VERSION 0)
set(Closure_MINOR_VERSION 0)
closure_print_status("Welcome! Closure version ${Closure_COLOR_BLUE}${Closure_WORLD_VERSION}.${Closure_MAJOR_VERSION}.${Closure_MINOR_VERSION}${Closure_COLOR_RESET}")

# common settings
set(Closure_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(Closure_LIB_PREFIX closure_)
set(Closure_CXX_FLAGS)

add_definitions(-DClosure_ROOT_DIR="${Closure_ROOT_DIR}")

closure_print_status("Closure_ROOT_DIR:${Closure_COLOR_BLUE}${Closure_ROOT_DIR}${Closure_COLOR_RESET}")
closure_print_status("Closure_LIB_PREFIX:${Closure_COLOR_BLUE}${Closure_LIB_PREFIX}${Closure_COLOR_RESET}")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(Closure_CXX_FLAGS " -stdlib=libstdc++ ")  # prevent default libc++ on new mac
endif()

## Build Options ===============================================================
option(Closure_WITH_UI "Build the library with UI" OFF)
option(Closure_WITH_TEST "Build Closure library with tests" ON)

# applications
option(Closure_WITH_APPS "Build the library with UI" OFF)

## General settings ============================================================

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

set(Closure_RUNTIME_OUTPUT_DIRECTORY ${Closure_ROOT_DIR}/bin CACHE PATH "Target for the binaries")
set(Closure_LIBRARY_OUTPUT_DIRECTORY ${Closure_ROOT_DIR}/lib CACHE PATH "Target for the libraries")
link_directories(${Closure_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${Closure_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${Closure_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Closure_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${Closure_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${Closure_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Closure_RUNTIME_OUTPUT_DIRECTORY})

# Detect OS and define macros appropriately
if(WIN32)
  add_definitions(-DWINDOWS)
  closure_print_status("Compiling for Windows")
elseif(ANDROID)
  add_definitions(-DANDROID)
  closure_print_status("Compiling for Android")
elseif(APPLE)
  add_definitions(-DAPPLE)
  MESSAGE(STATUS "Compiling for OSX")
elseif(UNIX)
  add_definitions(-DUNIX)
  MESSAGE(STATUS "Compiling for Unix")
else()
  closure_print_error("We don't support this platform!!!")
endif(WIN32)

set(CMAKE_CXX_FLAGS ${Closure_CXX_FLAGS})
closure_common_build_setting()

# Allow the developer to select if Dynamic or Static libraries are built
option(BUILD_SHARED_LIBS "Build Shared Libraries (preferred and required for the g2o plugin system)" OFF)
set(Closure_LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(Closure_LIB_TYPE SHARED)
  add_definitions(-DClosure_SHARED_LIBS)
endif()

## Dependencies configuration ==================================================
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/ClosureDependencies_Config.cmake")
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ClosureDependencies_Config.cmake.in"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
  file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/ClosureDependencies_Config.cmake.in"
              "${CMAKE_CURRENT_BINARY_DIR}/ClosureDependencies_Config.cmake")
endif()
include(ClosureDependencies_Config)
include(ClosureDependencies)

if(Closure_WITH_UI)
  if(NOT OPENGL_FOUND)
    closure_print_status("OpenGL not found! UI component and Applications require OpenGL")
  else()
    set(Closure_HAVE_OPENGL ON)
    closure_print_status("Found OpenGL!")
  endif()

  if(NOT GLUT_FOUND)
    closure_print_status("GLUT not found! UI component and Applications require GLUT")
  else()
    set(Closure_HAVE_GLUT ON)
    closure_print_status("Found GLUT!")
  endif()

  if(WIN32 AND NOT GLEW_FOUND)
    closure_print_status("GLEW not found! UI component and Applications require GLEW")
  elseif(GLEW_FOUND)
    set(Closure_HAVE_GLEW ON)
    closure_print_status("Found GLEW!")
  endif()
endif(Closure_WITH_UI)

## External libraries (included in this repo) ==================================
closure_print_status("CMAKE_CXX_FLAGS:${CMAKE_CXX_FLAGS}")
# Eigen
set(EIGEN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/closure_external)
add_subdirectory(external/eigen-3.2.1)
include_directories(external/eigen-3.2.1 external/eigen-3.2.1/unsupported)

# Look3D math
set(LOOK3D_LIB_PREFIX look3d_)
include_directories(external)
include_directories(external/look3d)
add_subdirectory(external/look3d/math)
list(APPEND Closure_EXTERNAL_LIBS look3d_math)

# openFabMap
add_subdirectory(external/openFabMap)

# FabMap
include_directories(external/openFabMap/include)
list(APPEND Closure_EXTERNAL_LIBS openFABMAP)

## Closure libraries ============================================================
set(Closure_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/closure
                         ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Closure_INCLUDE_DIRS})

set(Closure_LIBRARIES closure)
add_subdirectory(closure)  # where internal libraries are

## Apps ========================================================================
if(Closure_WITH_APPS)  # example applications
  add_subdirectory(apps)
endif(Closure_WITH_APPS)

# Installation =================================================================
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ClosureConfig.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/ClosureConfig.cmake"
               IMMEDIATE @ONLY)

# install FindClosure.cmake to common destination
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ClosureConfig.cmake
              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ClosureUtil.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)

# Uninstallation ===============================================================
# uninstall target
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
               IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)


# Useful paths configuration for windows 
if(WIN32)
  if(MSVC)
    if(CMAKE_CL_64)
      set(StructuralModeling_ARCH x64)
    else()
      set(StructuralModeling_ARCH x86)
    endif()
    if(MSVC_VERSION EQUAL 1400)
      set(StructuralModeling_RUNTIME vc8)
    elseif(MSVC_VERSION EQUAL 1500)
      set(StructuralModeling_RUNTIME vc9)
    elseif(MSVC_VERSION EQUAL 1600)
      set(StructuralModeling_RUNTIME vc10)
      get_filename_component(VC_IDE_PATH $ENV{VS100COMNTOOLS}/../IDE ABSOLUTE)
    elseif(MSVC_VERSION EQUAL 1700)
      set(StructuralModeling_RUNTIME vc11)
      get_filename_component(VC_IDE_PATH $ENV{VS110COMNTOOLS}/../IDE ABSOLUTE)
    elseif(MSVC_VERSION EQUAL 1800)
      set(StructuralModeling_RUNTIME vc12)
      get_filename_component(VC_IDE_PATH $ENV{VS120COMNTOOLS}/../IDE ABSOLUTE)
    endif()
  endif()
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/scripts/set_paths.bat.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/set_paths.bat")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/scripts/set_paths.bat.in" 
                 "${Closure_RUNTIME_OUTPUT_DIRECTORY}/set_paths.bat")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/scripts/set_paths_and_run_vc.bat.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/set_paths_and_run_vc.bat")
endif()
