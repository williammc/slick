cmake_minimum_required(VERSION 2.8)

project(slick_library_and_examples)

# Headers & sources ============================================================
set(sub_dirs geometry math scene util)
set(REG_EXT "[^/]*([.]cpp|[.]c|[.]cc|[.]h|[.]hpp)$")
foreach(sub_dir ${sub_dirs})
  file(GLOB slick_${sub_dir}_sources "${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir}/*.cc")
  file(GLOB slick_${sub_dir}_headers "${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir}/*.h")
  list(APPEND slick_sources ${slick_${sub_dir}_sources}
                            ${slick_${sub_dir}_headers})
  source_group("${sub_dir}" REGULAR_EXPRESSION "${CMAKE_CURRENT_SOURCE_DIR}/${sub_dir}/${REG_EXT}")
endforeach()

# Library ======================================================================
set(LIB_TARGET slick)
add_library(${LIB_TARGET} ${slick_sources})
set_property(TARGET ${LIB_TARGET} PROPERTY FOLDER "slick_libs")
set_property(TARGET ${LIB_TARGET} PROPERTY LINKER_LANGUAGE CXX)

#list(APPEND Slick_LIBS ${LIB_TARGET})

# Tests ========================================================================
if(Slick_WITH_TEST AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()  # has to be in root CMakeLists.txt for target 'test' to work!
  add_subdirectory(test)
endif(Slick_WITH_TEST AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)

## Examples ====================================================================
if(Slick_WITH_EXAMPLES)  # example applications
  add_subdirectory(examples)
endif(Slick_WITH_EXAMPLES)

