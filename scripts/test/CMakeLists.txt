# test/CMakeLists.txt

# Include the main source directory for header files
include_directories(${CMAKE_SOURCE_DIR}/src)
# file(GLOB eds_C "utils/*.cpp")
# file(GLOB eds_H "utils/*.h")
# file(GLOB eds ${eds_C} ${eds_H})

include_directories(~/include) #SDSL headers are here
LINK_DIRECTORIES(~/lib) #SDSL lib are here

# find_package(Boost REQUIRED)
# find_library(SDSL_LIB sdsl HINTS /home/draesdom/include/sdsl)
# Optionally, you can set the linker flags if needed

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -O3 -L/home/draesdom/lib -L/home/draesdom/Libs/boost_1_81_0/libs")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -g -L/home/draesdom/lib -L/home/draesdom/Libs/boost_1_81_0/libs")

# Gather all test files from the test directory
file(GLOB TEST_SOURCES "*.cpp")

# Create an executable for each test file and register it with CTest
foreach(TEST_SRC ${TEST_SOURCES})
    # Get the filename without the directory to use as the test name
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)

    # Create a test executable for each test source file
    add_executable(${TEST_NAME} ${TEST_SRC})    
    target_link_libraries(${TEST_NAME} "-lsdsl")    # Link with Boost libraries

    # Add the test executable to CTest
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()