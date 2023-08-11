cmake_minimum_required(VERSION 3.10)

# Project name and language
project(lineSurface LANGUAGES CXX)

# Set the build type to Debug or Release
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-g3 -Wall -Wextra -ansi -pedantic -Wno-unused-parameter -Wno-variadic-macros -Wconversion -std=c++0x)
    link_libraries(-g3)
else()
    add_compile_options(-g -O3)
    link_libraries(-g -O3)
endif()

# Include directories
include_directories(../../tools/)

# Find all .cpp files in the current directory
file(GLOB CPPFILES "*.cpp")

# Target name and source files
add_executable(lineSurface ${CPPFILES})

# If you have any libraries to link
# target_link_libraries(lineSurface ${LIBS})