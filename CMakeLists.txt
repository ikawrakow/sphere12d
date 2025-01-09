cmake_minimum_required(VERSION 3.14)
project("sphere12d" CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

if (NOT XCODE AND NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(circle circle.cpp randomGenerator.h randomGenerator.cpp)

add_executable(simple simple.cpp randomGenerator.h randomGenerator.cpp)

add_executable(smart smart.cpp randomGenerator.h randomGenerator.cpp)
