
[![Build Status](https://travis-ci.org/williammc/Slick.svg?branch=master)](https://travis-ci.org/williammc/Slick) : GCC4.8, Clang3.6 (Build & Tests)

Slick library and applications
Version 1.0.0

# About Slick 
`Slick` is a cool numerical library aiming to be concise, elegant, and fast.
The library includes only headers.
It is written in C++ and is complied with C++11 standard.
For testing, build with `Slick_WITH_TEST` flag on.

## Versions

### V1.0.0 (November 2014)
Initial code with Lie Group rigid transform (SO, SE), scene module (camera models), and geometry module (line, plane).

### [V1.1.0](https://github.com/williammc/Slick/releases/tag/v1.1) (April 2015)
Improved geometry module with more unittests for this module.


# Coding Guideline
`Slick` library and apps follows the C++ Style Guide from [googlecode.com](http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml)
For automatic code formating,  [clang-format](http://clang.llvm.org/docs/ClangFormat.html) can help.

# Dependencies (included in this repository)
- [Eigen-3.2+](www.eigen.tuxfamily.org) (included in ~/root/external): for matrix operation
- [gtest-1.7.0](code.google.com/p/googletest) (include in ~/root/external): for testing

# Testing
Build the tests by turn on `Slick_WITH_TEST` flag when doing cmake configuration.
To run all the test, run `ctest -C <config>`.
For instance, `ctest -C debug`

# Todo
- add more unittests (e.x.: for scene, util)
- more math :)