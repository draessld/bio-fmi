# src/CMakeLists.txt

# Find all .cpp files in src (excluding utils) and add them as executables
file(GLOB SRC_FILES CONFIGURE_DEPENDS "*.cpp")

# Collect all utility source files in src/utils
file(GLOB UTILS_FILES CONFIGURE_DEPENDS "utils/*.cpp")

# Define a utils library with the implementation files from src/utils
add_library(utils STATIC ${UTILS_FILES})

# Specify the include directory for utils headers
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utils)

# List of external libraries needed
set(EXTERNAL_LIBS sdsl boost_program_options divsufsort64 divsufsort)

foreach(SRC_FILE ${SRC_FILES})
    # Get the filename without directory path and extension
    get_filename_component(EXECUTABLE_NAME ${SRC_FILE} NAME_WE)

    # Add an executable for each .cpp file in src
    add_executable(${EXECUTABLE_NAME} ${SRC_FILE})

    # Link each executable to the utils library and external libraries
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE utils ${EXTERNAL_LIBS})
endforeach()

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# file(GLOB eds_C "utils/eds.cpp")
# file(GLOB eds_H "utils/eds.h")
# file(GLOB eds ${eds_C} ${eds_H})


# file(GLOB utils_H "utils/utils.hpp")


# Optionally, you can set the linker flags if needed


# file(GLOB SOURCES "*.cpp")

# # Create an executable for each .cpp file with a main function
# foreach(source_file ${SOURCES})
#     message(STATUS "Building executable: ${source_file}")
#     get_filename_component(file_name ${source_file} NAME_WE)
#     add_executable(${file_name} ${source_file} ${eds} ${utils_H})
#     target_link_libraries(${file_name} "-lsdsl")    # Link with Boost libraries
#     target_link_libraries(${file_name} "-lboost_program_options")    # Link with Boost libraries
#     target_link_libraries(${file_name} "-ldivsufsort64")    # Link with Boost libraries
#     target_link_libraries(${file_name} "-ldivsufsort")    # Link with Boost libraries
# endforeach()