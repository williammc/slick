include(ClosureUtil)
closure_print_status("Setup required dependencies for Closure tracking library")

##==============================================================================
## Boost
if(WIN32)
  set(Boost_USE_STATIC_LIBS  ON)
endif(WIN32)

find_package(Boost REQUIRED COMPONENTS system filesystem serialization)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(Closure_EXTERNAL_LIBS ${Closure_EXTERNAL_LIBS} ${Boost_LIBRARIES})
  link_directories(${Boost_LIBRARY_DIR})
  set(Closure_HAVE_SERIALIZATION ON)
  add_definitions(-DClosure_HAVE_SERIALIZATION)
  closure_print_status("Boost libs: ${Boost_LIBRARIES}")
endif(Boost_FOUND)

##==============================================================================
##  OpenCV
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
list(APPEND Closure_EXTERNAL_LIBS ${OpenCV_LIBRARIES})
closure_print_status("OpenCV include dirs:${OpenCV_INCLUDE_DIRS}")
closure_print_status("OpenCV libs:${OpenCV_LIBRARIES}")

##==============================================================================
## TBB
# set(TBB_ARCH_PLATFORM "ia32")
# set(TBB_COMPILER "vc12")
# set(TBB_INSTALL_DIR "D:/deps_msvc_common/tbb42_20140601oss")
find_package(TBB)
include_directories(${TBB_INCLUDE_DIR})
link_directories(${TBB_LIBRARY_DIRS})
list(APPEND EXTERNAL_LIBRARIES ${TBB_LIBRARIES})
set(TBB_BINARY_DIR ${TBB_INSTALL_DIR}/bin/${TBB_ARCH_PLATFORM}/${TBB_COMPILER})

closure_print_status("TBB_LIBRARIES:${TBB_LIBRARIES}")
closure_print_status("TBB_BINARY_DIR:${TBB_BINARY_DIR}")

if(TBB_FOUND)
  set(Closure_HAVE_TBB ON)
  add_definitions(-DClosure_HAVE_TBB)
endif()

##==============================================================================
# OpenGL
find_package(OpenGL)
closure_print_status("Found OpenGL ? ${OPENGL_FOUND}")
if(OPENGL_FOUND)
  closure_print_status("OpenGL INCLUDE: ${OPENGL_INCLUDE_DIR}")
  include_directories(${OPENGL_INCLUDE_DIR})
  list(APPEND Closure_EXTERNAL_LIBS ${OPENGL_LIBRARIES})
endif()

# GLUT
find_package(GLUT)
if(GLUT_FOUND)
  include_directories(${GLUT_INCLUDE_DIR})
  list(APPEND Closure_EXTERNAL_LIBS ${GLUT_LIBRARIES})
endif()

# GLEW
find_package(GLEW)
if(GLEW_FOUND)
  include_directories(${GLEW_INCLUDE_DIR})
  list(APPEND Closure_EXTERNAL_LIBS ${GLEW_LIBRARIES})
endif()