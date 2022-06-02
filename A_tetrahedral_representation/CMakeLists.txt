# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)


project(chanalyzer)


include_directories (BEFORE include)
include_directories(src)
include_directories(src/ConfigFile)
include_directories(src/jama)
include_directories(src/tnt)
include_directories(src/sturm)
include_directories(src/octree)


file(GLOB SRC
    "./src/*.h"
    "./src/*.cpp")


# CGAL and its components
find_package( CGAL )
include( ${CGAL_USE_FILE} )
add_definitions(-DENABLE_CGAL)

if (NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()


# Boost and its components
find_package( Boost REQUIRED COMPONENTS program_options system thread filesystem serialization REQUIRED )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )
include_directories (BEFORE include)
add_definitions(-DENABLE_BOOST_THREADS)

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()




add_executable(chanalyzer
							${SRC}
							./src/ConfigFile/ConfigFile.cpp
							./src/ConfigFile/ConfigFile.h 
							./src/sturm/sturm.cpp
							./src/sturm/util.cpp
							./src/sturm/solve.h
							./src/octree/octree.h
							./src/ply/ply.cpp
							./src/ply/ply.h 
							) 


# Creating entries for target: chanalyzer
# ############################

add_to_cached_list(CGAL_EXECUTABLE_TARGETS chanalyzer)

# Link the executable to CGAL and third-party libraries
target_link_libraries(chanalyzer PRIVATE CGAL::CGAL ${Boost_LIBRARIES})

# Set CMAKE_BUILD_TYPE to Release
set(CMAKE_BUILD_TYPE Release)
