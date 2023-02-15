cmake_minimum_required(VERSION 3.0)

project(algorithm-example)

# C++ Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(PROJECT_SOURCE_DIR source)

# Set Source Files
set(SOURCE ${PROJECT_SOURCE_DIR}/fast-multiply.cc)

add_executable(${PROJECT_NAME} ${SOURCE})