# test/CMakeLists.txt

# Enable testing
enable_testing()

# Specify the location of the utils headers and external libraries
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/utils)  # Include utils headers
include_directories(BEFORE ~/include)  # External headers
link_directories(BEFORE ~/lib)  # External libraries

# List of external libraries needed for tests
set(EXTERNAL_LIBS sdsl boost_program_options divsufsort64 divsufsort)

# Collect all test .cpp files in the test folder
file(GLOB TEST_FILES CONFIGURE_DEPENDS "*.cpp")

# Link the utils library and external libraries to each test
foreach(TEST_FILE ${TEST_FILES})
    # Get the test name from the filename without directory path and extension
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)

    # Add an executable for each test file
    add_executable(${TEST_NAME} ${TEST_FILE})

    # Link each test executable to the utils library and external libraries
    target_link_libraries(${TEST_NAME} PRIVATE utils ${EXTERNAL_LIBS})

    # Add each test to CTest
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
