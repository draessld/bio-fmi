# For tests

macro (create_test TESTNAME)
    add_executable (${TESTNAME} ${TESTNAME}.cpp)
    target_include_directories (${TESTNAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries (${TESTNAME} gtest_main)
    add_test (NAME ${TESTNAME} COMMAND ${TESTNAME})
endmacro (create_test)
