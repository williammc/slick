include(cotire)
## Nice colors in command ======================================================
# on Windows, use ConEmu console emulator to have ansi colors output
# recommended to use with Clink (for useful commands)
string(ASCII 27 Esc)
set(CMAKER_COLOR_NORMAL          "")
set(CMAKER_COLOR_RESET           "${Esc}[m")
set(CMAKER_COLOR_BOLD            "${Esc}[1m")
set(CMAKER_COLOR_RED             "${Esc}[31m")
set(CMAKER_COLOR_GREEN           "${Esc}[32m")
set(CMAKER_COLOR_YELLOW          "${Esc}[33m")
set(CMAKER_COLOR_BLUE            "${Esc}[34m")
set(CMAKER_COLOR_MAGENTA         "${Esc}[35m")
set(CMAKER_COLOR_CYAN            "${Esc}[36m")
set(CMAKER_COLOR_BOLD_RED        "${Esc}[1;31m")
set(CMAKER_COLOR_BOLD_GREEN      "${Esc}[1;32m")
set(CMAKER_COLOR_BOLD_YELLOW     "${Esc}[1;33m")
set(CMAKER_COLOR_BOLD_BLUE       "${Esc}[1;34m")
set(CMAKER_COLOR_BOLD_MAGENTA    "${Esc}[1;35m")
set(CMAKER_COLOR_BOLD_CYAN       "${Esc}[1;36m")
set(CMAKER_COLOR_BG_RED          "${Esc}[41m")
set(CMAKER_COLOR_BG_GREEN        "${Esc}[42m")
set(CMAKER_COLOR_BG_YELLOW       "${Esc}[43m")
set(CMAKER_COLOR_BG_BLUE         "${Esc}[44m")
set(CMAKER_COLOR_BG_MAGENTA      "${Esc}[45m")
set(CMAKER_COLOR_BG_CYAN         "${Esc}[46m")

## Useful commands  ============================================================
# Cull out paths in given list of libs and keep names
macro(cmaker_cull_library_paths LIBS)
  set(IN_LIBS ${${LIBS}})
  set(${LIBS} "")
  foreach(CULL_LIB ${IN_LIBS})
    get_filename_component(FN "${CULL_LIB}" NAME)
    list(APPEND ${LIBS} ${FN})
  endforeach()
endmacro(cmaker_cull_library_paths)

# To add apps folder from external directories
macro(cmaker_add_external_apps_path in_directory in_binary_path in_message)
  if(EXISTS ${in_directory})
    message(STATUS ${in_message})
    add_subdirectory(${in_directory}  ${in_binary_path})
  else()
    message(WARNING "Input directory ${in_directory} not exists")
  endif()
endmacro(cmaker_add_external_apps_path)

# To print nice message with information about the source file
# Output: file_name:: message
macro(cmaker_print_status_colorized color msg)
  get_filename_component(CURRENT_FILENAME "${CMAKE_CURRENT_LIST_FILE}" NAME)
  get_filename_component(CURRENT_PATH "${CMAKE_CURRENT_LIST_FILE}" PATH)
  get_filename_component(CURRENT_PARENT_DIRECTORY "${CURRENT_PATH}" NAME)
  message(STATUS "${color}${CURRENT_PARENT_DIRECTORY}/${CURRENT_FILENAME}:${CMAKE_CURRENT_LIST_LINE} >> ${msg}${CMAKER_COLOR_RESET}")
endmacro(cmaker_print_status_colorized)

# To print nice message with information about the source file
# Output: file_name:: message
macro(cmaker_print_status msg)
  get_filename_component(CURRENT_FILENAME "${CMAKE_CURRENT_LIST_FILE}" NAME)
  get_filename_component(CURRENT_PATH "${CMAKE_CURRENT_LIST_FILE}" PATH)
  get_filename_component(CURRENT_PARENT_DIRECTORY "${CURRENT_PATH}" NAME)
  if(NOT ARGN)
    message(STATUS "${CMAKER_COLOR_GREEN}${CURRENT_PARENT_DIRECTORY}/${CURRENT_FILENAME}:${CMAKE_CURRENT_LIST_LINE} >> ${msg}${ARGN}${CMAKER_COLOR_RESET}")
  else()
    message(STATUS "${CMAKER_COLOR_GREEN}${CURRENT_PARENT_DIRECTORY}/${CURRENT_FILENAME}:${CMAKE_CURRENT_LIST_LINE} >> \n${CMAKER_COLOR_RESET}")
    list(APPEND msg ${ARGN})
    foreach(m in ${msg})
      message(STATUS "${CMAKER_COLOR_GREEN}${m}${CMAKER_COLOR_RESET}")
    endforeach()
endif()
endmacro(cmaker_print_status)

# To print nice error message with information about the source file
macro(cmaker_print_error msg) # TODO: improve this for ARGN
  get_filename_component(CURRENT_FILENAME "${CMAKE_CURRENT_LIST_FILE}" NAME)
  get_filename_component(CURRENT_PATH "${CMAKE_CURRENT_LIST_FILE}" PATH)
  get_filename_component(CURRENT_PARENT_DIRECTORY "${CURRENT_PATH}" NAME)
  message(FATAL_ERROR "${CMAKER_COLOR_RED}${CURRENT_PARENT_DIRECTORY}/${CURRENT_FILENAME}:${CMAKE_CURRENT_LIST_LINE} >> ${msg}${ARGN}${CMAKER_COLOR_RESET}")
endmacro(cmaker_print_error)

# Contains common (useful) cmake macros or routines ============================
# useful routine for building a binary
# requires ${ALL_LIBRARRIES}

macro(cmaker_add_library target lib_type)
  add_library(${target} ${lib_type} ${ARGN})
  if(ALL_LIBRARIES)
    add_dependencies(${target} ${ALL_LIBRARIES})
    target_link_libraries(${target} ${ALL_LIBRARIES})
  else()
    cmaker_print_status("ALL_LIBRARIES is not set for target:${target}")
  endif()
endmacro()

macro(cmaker_add_library_cotire target lib_type)
  add_library(${target} ${lib_type} ${ARGN})
  if(ALL_LIBRARIES)
    add_dependencies(${target} ${ALL_LIBRARIES})
    target_link_libraries(${target} ${ALL_LIBRARIES})
    set_target_properties(${target} PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY")
    cotire(${target})
  else()
    cmaker_print_status("ALL_LIBRARIES is not set for target:${target}")
  endif()
endmacro()

macro(cmaker_add_executable target)
  set(srcs ${target}.cc ${ARGN})
  add_executable(${target} ${srcs})
  if(ALL_LIBRARIES)
    add_dependencies(${target} ${ALL_LIBRARIES})
    target_link_libraries(${target} ${ALL_LIBRARIES})
  else()
    cmaker_print_status("ALL_LIBRARIES is not set for target:${target}")
  endif()
endmacro(cmaker_add_executable)

macro(cmaker_add_executable_cotire target)
  set(srcs ${target}.cc ${ARGN})
  add_executable(${target} ${srcs})
  if(ALL_LIBRARIES)
    add_dependencies(${target} ${ALL_LIBRARIES})
    target_link_libraries(${target} ${ALL_LIBRARIES})
    set_target_properties(${target} PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY")
    cotire(${target})
  else()
    cmaker_print_status("ALL_LIBRARIES is not set for target:${target}")
  endif()
endmacro(cmaker_add_executable_cotire)

# add a test
macro(cmaker_add_test target)
  cmaker_add_executable(${target} ${ARGN})
  add_test(NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro(cmaker_add_test)

macro(cmaker_add_test_cotire target)
  cmaker_add_executable(${target} ${ARGN})
  add_test(NAME ${target} COMMAND $<TARGET_FILE:${target}>)
  set_target_properties(${target} PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY")
  cotire(${target})
endmacro(cmaker_add_test_cotire)

# commong build settings (mostly to deal with Windows )
macro(cmaker_common_build_setting)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
    option(ENABLE_NEON "Enable NEON" ON)
    if(ENABLE_NEON)
      add_definitions(-DENABLE_NEON)  # enable NEON code for ARM
    endif()
  else()
    option(ENABLE_SSE "Enable SSE" ON)
    if(ENABLE_SSE)
      add_definitions(-DENABLE_SSE)   # SSE code
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${SSE_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${SSE_FLAGS}")
  endif()
  
  add_definitions(-DM_SQRT1_2=0.707106781186547524401)  # missing in MSVC2013
    
  if(WIN32 AND MSVC)
    #add_definitions(/DEIGEN_DONT_ALIGN_STATICALLY)  # this is unusual
    add_definitions(/D_CRT_SECURE_NO_WARNINGS)
    add_definitions(/DWIN32_LEAN_AND_MEAN)  # speed up compiler by omitting some headers
    add_definitions(/DM_SQRT1_2=0.707106781186547524401)  # missing in MSVC2013
    add_definitions(/DM_PI=3.14159265358979323846) # missing in MSVC2013
    add_definitions("/MP")  # multiple processes compilation
    
    add_definitions(/DNOMINMAX)  # resolve conflicts of std::min()/std::max() on Windows MSVC
    add_definitions(/D_USE_MATH_DEFINES)
    if(BUILD_SHARED_LIBS)  # disable warning on missing DLL interfaces
      add_definitions("/wd4251")
    endif()

    #add_definitions(/D_HAS_EXCEPTIONS=0)  # disable exception handling code
    #add_definitions(/D"__uncaught_exception()=true")
    #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NOENTRY ")  # this causes IsValidPtr errors
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qpar " )  # parallel code generation
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /D_SECURE_SCL=0 /GR- /fp:fast  /GS- /W0 ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /D_SECURE_SCL=0 /fp:fast  /GS- /W0 ")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /D_SECURE_SCL=0 ")

  else(WIN32 AND MSVC)
    if(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
      # SSE is not supported on ARM!
      if(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        add_definitions(-mfpu=neon -mfloat-abi=softfp -march=armv7-a)  # vectorization for ARM
      else()
        #add_definitions(-march=native)
      endif()
    endif()
  endif(WIN32 AND MSVC)
endmacro(cmaker_common_build_setting)

# Detect OS and define macros appropriately
macro(cmaker_detect_os)
if(WIN32)
  add_definitions(/DWINDOWS /D_WIN32)
  cmaker_print_status("Compiling for Windows")
elseif(ANDROID)
  add_definitions(-DANDROID)
  cmaker_print_status("Compiling for Android")
elseif(APPLE)
  add_definitions(-DAPPLE)
  MESSAGE(STATUS "Compiling for OSX")
elseif(UNIX)
  add_definitions(-DUNIX)
  MESSAGE(STATUS "Compiling for Unix")
else()
  cmaker_print_error("We don't support this platform!!!")
endif(WIN32)
endmacro(cmaker_detect_os)
