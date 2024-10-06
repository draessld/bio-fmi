# src/CMakeLists.txt

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

file(GLOB eds_C "utils/eds.cpp")
file(GLOB eds_H "utils/eds.h")
file(GLOB eds ${eds_C} ${eds_H})


file(GLOB utils_H "utils/utils.hpp")

include_directories(~/include) #SDSL headers are here
LINK_DIRECTORIES(~/lib) #SDSL lib are here

# Optionally, you can set the linker flags if needed

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -O3 -L/home/draesdom/lib -L/home/draesdom/Libs/boost_1_81_0/libs")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -g -L/home/draesdom/lib -L/home/draesdom/Libs/boost_1_81_0/libs")

file(GLOB SOURCES "*.cpp")

# Create an executable for each .cpp file with a main function
foreach(source_file ${SOURCES})
    message(STATUS "Building executable: ${source_file}")
    get_filename_component(file_name ${source_file} NAME_WE)
    add_executable(${file_name} ${source_file} ${eds} ${utils_H})
    target_link_libraries(${file_name} "-lsdsl")    # Link with Boost libraries
    target_link_libraries(${file_name} "-lboost_program_options")    # Link with Boost libraries
    target_link_libraries(${file_name} "-ldivsufsort64")    # Link with Boost libraries
    target_link_libraries(${file_name} "-ldivsufsort")    # Link with Boost libraries
endforeach()