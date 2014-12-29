cmake_minimum_required(VERSION 2.8.9)

set(THE_PROJECT_NAME Slick)
project(${THE_PROJECT_NAME})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)  # for organizing code to folders

# Look for our CMake modules paths
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakerCommonUtil.cmake)
  set(CMAKER_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
else()
find_path(CMAKER_ROOT NAMES cmake/CMakerCommonUtil.cmake
            PATHS ${CMAKE_CURRENT_SOURCE_DIR}
                  $ENV{CMAKER_ROOT})
endif()
if(NOT CMAKER_ROOT)
  message(ERROR "Environment variable CMAKER_ROOT not set!")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKER_ROOT}/cmake
                              ${CMAKER_ROOT}/cmake/modules
                              ${CMAKE_CURRENT_SOURCE_DIR}/cmake
                              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules
                              ${CMAKE_CURRENT_BINARY_DIR})

set(Slick_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
if(NOT THE_PROJECT_ROOT OR THE_PROJECT_ROOT STREQUAL "")
  set(THE_PROJECT_ROOT ${Slick_ROOT})
endif()

include(CMakerCommonUtil)  # include useful macros
include(CMakerCommonSettings)

set(Slick_WORLD_VERSION 1)
set(Slick_MAJOR_VERSION 0)
set(Slick_MINOR_VERSION 0)
cmaker_print_status("Welcome! Slick version ${Slick_COLOR_BLUE}${Slick_WORLD_VERSION}.${Slick_MAJOR_VERSION}.${Slick_MINOR_VERSION}${Slick_COLOR_RESET}")
cmaker_print_status("Slick cmake module paths: ${Slick_COLOR_BLUE}${CMAKE_MODULE_PATH}${Slick_COLOR_RESET}")
cmaker_detect_os()

set(Slick_CXX_FLAGS)
add_definitions(-DSlick_ROOT="${Slick_ROOT}")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(Slick_CXX_FLAGS " -stdlib=libstdc++ ")  # prevent default libc++ on new mac
endif()

## Build Options ===============================================================
option(Slick_WITH_TEST "Build Slick library with tests" OFF)
option(Slick_WITH_EXAMPLES "Build the library with UI" OFF)

set(CMAKE_CXX_FLAGS ${Slick_CXX_FLAGS})
cmaker_common_build_setting()

## External libraries (included in this repo) ==================================
# Eigen
set(EIGEN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/slick_external)
#add_subdirectory(external/eigen-3.2.2)
include_directories(external/eigen-3.2.2 external/eigen-3.2.2/unsupported)

# GTest
if(Slick_WITH_TEST AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/gtest-1.7.0)
  add_subdirectory(external/gtest-1.7.0)
  include_directories(external/gtest-1.7.0/include)
  set_property(TARGET gtest PROPERTY FOLDER "external/gtest")
  set_property(TARGET gtest_main PROPERTY FOLDER "external/gtest")
endif()

## Slick library ===============================================================
set(Slick_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Slick_INCLUDE_DIRS})

add_subdirectory(slick)  # where internal libraries are
