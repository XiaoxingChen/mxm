cmake_minimum_required(VERSION 2.14)
project(mxm)


if(CMAKE_SYSTEM_NAME MATCHES "^Android")
    add_compile_options(-std=c++17)
    set(THREAD_LIB stdc++)
elseif(CMAKE_SYSTEM_NAME MATCHES "^Linux")
    add_compile_options(-std=c++17)
    set(THREAD_LIB pthread)
elseif(CMAKE_SYSTEM_NAME MATCHES "^Windows")
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
    endif(MSVC)
else()
    add_compile_options(-std=c++17)
endif()

include_directories(inc)

# add_library(mxm STATIC src/mxm.cpp)
add_library(mxm STATIC src/mxm.cpp)
target_compile_definitions(mxm PUBLIC MXM_COMPILED_LIB)

if(BUILD_TEST)
    add_subdirectory(tests)
endif(BUILD_TEST)
