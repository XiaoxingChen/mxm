cmake_minimum_required(VERSION 3.2)

set(TEST_MAIN_SRC
    test_main.cpp
    test_linalg.cpp
    test_cv_basic.cpp
    test_lie_alg.cpp
    test_joint.cpp)

add_executable(test_main ${TEST_MAIN_SRC})
target_include_directories(test_main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
target_compile_definitions(test_main PUBLIC MXM_COMPILED_LIB)
target_link_libraries(test_main PUBLIC mxm)
