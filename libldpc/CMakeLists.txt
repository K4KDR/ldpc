cmake_minimum_required(VERSION 3.0)

##############################################
# Create target and set properties

set(ldpc_SOURCES
    src/decoder.cpp
    src/encoder.cpp
    src/ldpc.cpp
)

add_library(ldpc SHARED ${ldpc_SOURCES})

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(ldpc::ldpc ALIAS ldpc)

target_include_directories(ldpc
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Set compiler to strict mode when compiling this binary
target_compile_options(ldpc
    PRIVATE
        $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
            -Werror -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion>
        $<$<C_COMPILER_ID:MSVC>:
            /WX /W4>
)

##############################################
# Installation instructions

install(DIRECTORY include/ DESTINATION include)

install(TARGETS ldpc
    EXPORT ldpcTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)
