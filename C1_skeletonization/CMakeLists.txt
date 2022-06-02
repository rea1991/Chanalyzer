cmake_minimum_required(VERSION 3.1...3.15)

project(skeletonization)

include_directories (BEFORE include)
include_directories(src)
include_directories(src/auxiliary)

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

add_executable(skeletonization
							${SRC}
) 

# Link the executable to CGAL and third-party libraries
target_link_libraries(skeletonization PRIVATE CGAL::CGAL ${Boost_LIBRARIES})

# Set CMAKE_BUILD_TYPE to Release
set(CMAKE_BUILD_TYPE Release)
