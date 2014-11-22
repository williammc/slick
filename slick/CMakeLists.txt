cmake_minimum_required(VERSION 2.8)

project(look3d_math)

# Headers & sources ============================================================
file(GLOB look3d_math_sources "*.cc")
file(GLOB look3d_math_headers "*.h")

# Library ======================================================================
set(LIB_TARGET ${LOOK3D_LIB_PREFIX}math)
set(LIB_DEPENDENCIES ${LOOK3D_LIB_PREFIX}core)
add_library(${LIB_TARGET} ${LOOK3D_LIB_TYPE} ${look3d_math_sources}
                                             ${look3d_math_headers})
target_link_libraries(${LIB_TARGET} ${LIB_DEPENDENCIES} ${LOOK3D_EXTERNAL_LIBS})
set_property(TARGET ${LIB_TARGET} PROPERTY FOLDER "look3d_libs")

list(APPEND LOOK3D_LIBRARIESS ${LIB_TARGET})

# Tests ========================================================================
if(LOOK3D_WITH_TEST)
  enable_testing()  # has to be in root CMakeLists.txt for target 'test' to work!
  add_subdirectory(test)
endif(LOOK3D_WITH_TEST)

# Installation =================================================================
install(TARGETS ${LIB_TARGET}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${look3d_math_headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/look3d/math)
