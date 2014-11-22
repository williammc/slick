// Description: import/export macros for creating DLLS with MSVC
// Use SM_<module>_API macro to declare it.
// Sparate macros for each DLL required.
#ifndef LOOK3D_MATH_API_H
#define LOOK3D_MATH_API_H

#ifdef _MSC_VER  // Microsoft compiler:
  #ifdef LOOK3D_SHARED_LIBS
    #ifdef look3d_math_EXPORTS
      #define LOOK3D_MATH_API __declspec(dllexport)
    #else
      #define LOOK3D_MATH_API __declspec(dllimport)
    #endif
  #else
    #define LOOK3D_MATH_API
  #endif

#else  // not MSVC
  #define LOOK3D_MATH_API
#endif

#endif  // LOOK3D_MATH_API_H
