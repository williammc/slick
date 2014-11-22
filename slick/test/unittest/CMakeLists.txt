cmake_minimum_required(VERSION 2.8)
project(slick_unit_testing)

# common settings
cmaker_common_build_setting()

##==============================================================================
#find_package(GTest)
#cmaker_print_status("GTEST_FOUND:${GTEST_FOUND}")

# cmaker_print_status("GTEST_INCLUDE_DIRS:${GTEST_INCLUDE_DIRS}")
# if(NOT GTEST_FOUND)
#   return()
# endif()

# include_directories(${GTEST_INCLUDE_DIRS})
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

set(ALL_LIBRARIES ${Slick_LIBS} ${GTEST_BOTH_LIBRARIES})

## Tests =======================================================================
file(GLOB test_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
foreach(filename ${test_sources})
  get_filename_component(test_name ${filename} NAME_WE)
  cmaker_add_test(${test_name} ${test_name}.cc util.h)
  set_property(TARGET ${test_name} PROPERTY FOLDER "slick_libs/tests/unittest")
endforeach()