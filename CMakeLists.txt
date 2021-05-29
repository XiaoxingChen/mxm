cmake_minimum_required(VERSION 2.14)
project(mxm)
add_compile_options(-std=c++14)
set(CMAKE_BUILD_TYPE Debug)

if(CMAKE_SYSTEM_NAME MATCHES "^Android")
    set(THREAD_LIB stdc++)
elseif(CMAKE_SYSTEM_NAME MATCHES "^Linux")
    set(THREAD_LIB pthread)
endif()

include_directories(inc)

# add_library(mxm STATIC src/mxm.cpp)
add_library(mxm STATIC src/mxm.cpp)
target_compile_definitions(mxm PUBLIC MXM_COMPILED_LIB)

if(BUILD_TEST)
    add_subdirectory(tests)
endif(BUILD_TEST)