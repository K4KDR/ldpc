############################################################
# Compute generator matrix
############################################################

option (LDPC_BUILD_GENERATOR "Build the ldpc_compute_generator binary (requires NTL)" OFF)

if (LDPC_BUILD_GENERATOR)
	message(STATUS "Optional ldpc_compute_generator binary will be build. Disable this with -DLIBLDPC_BUILD_GENERATOR=Off")

	find_package(NTL REQUIRED MODULE)
	
	add_executable(ldpc_compute_generator ldpc_compute_generator.cpp)
	target_link_libraries(ldpc_compute_generator NTL::NTL)
	install(TARGETS ldpc_compute_generator DESTINATION bin)
	
	# Set compiler to strict mode when compiling this binary
    target_compile_options(ldpc_compute_generator
        PRIVATE
            $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
                -Werror -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion>
            $<$<C_COMPILER_ID:MSVC>:
                /WX /W4>
    )
else()
	message(STATUS "Optional ldpc_compute_generator binary is not build. Enable this with -DLIBLDPC_BUILD_GENERATOR=On")
endif()
