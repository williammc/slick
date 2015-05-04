project(slick_test)

# common settings
cmaker_common_build_setting()

##==============================================================================
set(ALL_LIBRARIES ${Slick_LIBS} ${Slick_EXTERNAL_LIBS})
##==============================================================================
##  executables
file(GLOB test_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
foreach(filename ${test_sources})
  get_filename_component(test_name ${filename} NAME_WE)
  cmaker_add_test(${test_name} ${test_name}.cc)
  set_property(TARGET ${test_name} PROPERTY FOLDER "slick_tests")
endforeach()

add_subdirectory(unittest)
