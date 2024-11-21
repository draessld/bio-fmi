# Root CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Set Debug-specific flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -g -lboost_program_options")

# Set default build type to Debug if not specified
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)



include_directories(BEFORE ~/include)
link_directories(BEFORE ~/lib)

enable_testing()

# Add the src directory to the project
add_subdirectory(src)
add_subdirectory(test)
