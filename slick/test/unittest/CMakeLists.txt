cmake_minimum_required(VERSION 2.8)
project(Look3D_math_unit_testing)

# common settings
look3d_common_build_setting()

##==============================================================================
find_package(GTest)
look3d_print_status("GTEST_FOUND:${GTEST_FOUND}")

look3d_print_status("GTEST_INCLUDE_DIRS:${GTEST_INCLUDE_DIRS}")
if(NOT GTEST_FOUND)
  return()
endif()

include_directories(${GTEST_INCLUDE_DIRS})

set(ALL_LIBRARIES ${Look3D_LIBRARIES} ${LOOK3D_EXTERNAL_LIBS}
                  ${GTEST_BOTH_LIBRARIES})

## Tests =======================================================================
set(math_tests so2_unittest so3_unittest se2_unittest se3_unittest)
foreach(test ${math_tests})
  look3d_add_test(${test})
  set_property(TARGET ${test} PROPERTY FOLDER "look3d_tests/math/unittest")
endforeach()

