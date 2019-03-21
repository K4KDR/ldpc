cmake_minimum_required(VERSION 3.0)

project(ldpc)

add_subdirectory(libldpc)
add_subdirectory(bin)

enable_testing()
add_subdirectory(tests)
