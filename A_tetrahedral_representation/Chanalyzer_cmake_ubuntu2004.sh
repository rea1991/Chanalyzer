#!/bin/bash

#cgal_create_CMakeLists -s main_test
cd ./build
cmake .. #-DCMAKE_CXX_COMPILER=g++ #-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
make
