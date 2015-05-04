cmake_minimum_required(VERSION 2.8.7)

set(THE_PROJECT_NAME Slick)
project(${THE_PROJECT_NAME})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)  # for organizing code to folders

# Look for our CMake modules paths
set(CMAKER_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/cmaker)
if(NOT EXISTS "${CMAKER_ROOT}/cmake/CMakerCommonUtil.cmake")
  message(FATAL_ERROR "cmaker folder not exist")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKER_ROOT}/cmake
                              ${CMAKER_ROOT}/cmake/modules
                              ${CMAKE_CURRENT_SOURCE_DIR}/cmake
                              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules
                              ${CMAKE_CURRENT_BINARY_DIR})

set(Slick_ROOT ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
if(NOT THE_PROJECT_ROOT OR THE_PROJECT_ROOT STREQUAL "")
  set(THE_PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
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
  set(Slick_CXX_FLAGS " -stdlib=libc++ ")  # prevent default libc++ on new mac
endif()

if(NOT WIN32)
  set(Slick_CXX_FLAGS " ${Slick_CXX_FLAGS} -std=c++11 ") 
endif()
cmaker_print_status("Slick_CXX_FLAGS:${Slick_CXX_FLAGS}")

## Build Options ===============================================================
option(Slick_WITH_TEST "Build Slick library with tests" OFF)
if(Slick_WITH_TEST_ON)
  set(Slick_WITH_TEST ON)
endif()
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
  enable_testing()
  add_subdirectory(external/gtest-1.7.0)
  include_directories(external/gtest-1.7.0/include)
  set_property(TARGET gtest PROPERTY FOLDER "external/gtest")
  set_property(TARGET gtest_main PROPERTY FOLDER "external/gtest")
endif()

## Slick library ===============================================================
set(Slick_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Slick_INCLUDE_DIRS})

add_subdirectory(slick)  # where internal libraries are


# useful configuration (e.x.: put scripts in right places) =====================
if(WIN32)
  configure_file("${CMAKER_ROOT}/scripts/run_ninja_build.bat.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/run_ninja_build.bat")
  if(NOT EXISTS "${THE_LIB_RUNTIME_OUTPUT_DIRECTORY}/run_windbg.bat")
    configure_file("${CMAKER_ROOT}/scripts/run_windbg.bat.in" 
                 "${THE_LIB_RUNTIME_OUTPUT_DIRECTORY}/run_windbg.bat")
  endif()
  set(THE_BUILD_CMD "run_ninja_build.bat")
  set(THE_FILE_REGEX "^(.*)\\\\(([0-9]*)\\\\)(\\\\s:\\\\s).*error(.*)")
else()
  # configure it to a temp location
  configure_file("${CMAKER_ROOT}/scripts/run_make_build.sh.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/run_make_build.sh")
  # now copy the temporary into the final destination, setting the permissions
  file(COPY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/run_make_build.sh"
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
          FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
          GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  set(THE_BUILD_CMD "run_make_build.sh")
  set(THE_FILE_REGEX "^(.*):([0-9]+):([0-9]+).*error(.*)")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/scripts/Slick.sublime-project"
               "${CMAKE_CURRENT_BINARY_DIR}/Slick.sublime-project" @ONLY)

