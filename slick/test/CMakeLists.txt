cmake_minimum_required(VERSION 2.8)
project(look3d_math_test)

# common settings
look3d_common_build_setting()

##==============================================================================
set(ALL_LIBRARIES ${LOOK3D_LIBRARIES} ${LOOK3D_EXTERNAL_LIBS})
##==============================================================================
##  executables

set(math_tests test_se2 test_se3 test_sl test_so2 test_so3)
foreach(test ${math_tests})
  look3d_add_test(${test} ${test}.cc)
  set_property(TARGET ${test} PROPERTY FOLDER "look3d_tests/math")
endforeach()

add_subdirectory(unittest)
