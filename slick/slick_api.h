// Description: import/export macros for creating DLLS with MSVC
// Use SLICK_API macro to declare it.
// Sparate macros for each DLL required.
#pragma once

#ifdef _MSC_VER  // Microsoft compiler:
  #ifdef Slick_SHARED_LIBS
    #ifdef slick_EXPORTS
      #define SLICK_API __declspec(dllexport)
    #else
      #define SLICK_API __declspec(dllimport)
    #endif
  #else
    #define SLICK_API
  #endif

#else  // not MSVC
  #define SLICK_API
#endif

