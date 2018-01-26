#[[ Set CMake minimum required version ]]
cmake_minimum_required(VERSION 3.6)

#[[ Set C++ standard version ]]
set(CMAKE_CXX_STANDARD 14)

#[[ Set header files ]]
set(HEADER_FILES
        plane.h)

#[[ Set source files ]]
set(SOURCE_FILES
        plane.cpp)

#[[ Set test entry point ]]
set(ENTRY_POINT
        main.cpp)

#[[ Add executable target ]]
add_executable(plane-example
        ${HEADER_FILES}
        ${SOURCE_FILES}
        ${ENTRY_POINT})

#[[ Add include directories ]]
target_include_directories(plane-example PUBLIC ../../Source)
