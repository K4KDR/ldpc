#############################################
# Create target and set properties

set(ldpc_SOURCES
    include/ldpc/decoder.h
    include/ldpc/encoder.h
    include/ldpc/ldpc.h
    src/decoder.cpp
    src/encoder.cpp
    src/ldpc.cpp
)

add_library(ldpc SHARED ${ldpc_SOURCES})

target_include_directories(ldpc
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}/
)

# Make all symbols hidden by default and explicitly export them afterwards
set_property(TARGET ldpc PROPERTY C_VISIBILITY_PRESET hidden)
set_property(TARGET ldpc PROPERTY VISIBILITY_INLINES_HIDDEN 1)


include(GenerateExportHeader)

generate_export_header(ldpc
    EXPORT_FILE_NAME ldpc/ldpc_export.h
    EXPORT_MACRO_NAME LDPC_EXPORT
    DEPRECATED_MACRO_NAME LDPC_DEPRECATED
    NO_EXPORT_MACRO_NAME LDPC_NO_EXPORT
    NO_DEPRECATED_MACRO_NAME LDPC_NO_EXPORT_DEPRECATED
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

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ldpc/ldpc_export.h DESTINATION include/ldpc
)
install(TARGETS ldpc
    EXPORT ldpcTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(ldpc::ldpc ALIAS ldpc)
