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
