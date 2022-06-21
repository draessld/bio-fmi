cmake_minimum_required(VERSION 3.0.0)
project(bio-fmi VERSION 0.1.0)

include(CTest)
enable_testing()


set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../exe/")

# ############################################################
set(BIO-FMI-BUILD "bio-fmi-build")
set(BIO-FMI-LOCATE "bio-fmi-locate")

set(BIO-FMI-BUILD_EXEC "src/bio-fmi-build.cpp")
set(BIO-FMI-LOCATE_EXEC "src/bio-fmi-locate.cpp")

# GLOBAL OPTIONS
############################################################
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -O3 -lboost_program_options")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++17 -g -lboost_program_options")

# files for BIO-FMI
############################################################
file(GLOB BIO-FMI_C "src/index/*.cpp")
file(GLOB BIO-FMI_H "src/index/*.h")
file(GLOB BIO-FMI ${BIO-FMI_C} ${BIO-FMI_H})

# # BIO-FMI UTILS
# ############################################################
file(GLOB_RECURSE UTILS "src/utils/*.cpp")

# SDSL
############################################################
include_directories(~/include) #SDSL headers are here
LINK_DIRECTORIES(~/lib) #SDSL lib are here

# EXECUTABLES
############################################################

# bio-fmi
##############################
add_executable(${BIO-FMI-BUILD} ${BIO-FMI-BUILD_EXEC} ${UTILS} ${BIO-FMI})
target_link_libraries(${BIO-FMI-BUILD} "-lsdsl")
target_link_libraries(${BIO-FMI-BUILD} "-ldivsufsort")
target_link_libraries(${BIO-FMI-BUILD} "-ldivsufsort64")
target_link_libraries(${BIO-FMI-BUILD} "-lboost_program_options")

add_executable(${BIO-FMI-LOCATE} ${BIO-FMI-LOCATE_EXEC} ${UTILS} ${BIO-FMI})
target_link_libraries(${BIO-FMI-LOCATE} "-lsdsl")
target_link_libraries(${BIO-FMI-LOCATE} "-ldivsufsort")
target_link_libraries(${BIO-FMI-LOCATE} "-ldivsufsort64")
target_link_libraries(${BIO-FMI-LOCATE} "-lboost_program_options")
