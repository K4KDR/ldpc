cmake_minimum_required(VERSION 3.0)

project(ldpc)

add_subdirectory(libldpc)
add_subdirectory(bin)

option(LIBLDPC_UNITTESTS "Enable unittests" Off)

if(LIBLDPC_UNITTESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
