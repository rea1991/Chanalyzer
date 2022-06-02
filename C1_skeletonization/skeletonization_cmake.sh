#!/bin/bash

#cgal_create_CMakeLists -s ./src/main
cd ./build
cmake .. -DCMAKE_OSX_ARCHITECTURES=x86_64 -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++
make
