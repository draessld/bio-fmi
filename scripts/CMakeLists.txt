# Root CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

# Add the src directory to the project
add_subdirectory(src)
add_subdirectory(test)
