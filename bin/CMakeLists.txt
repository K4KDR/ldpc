cmake_minimum_required(VERSION 3.1)

project(ldpc-bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/bin/cmake/Modules/")

############################################################
# Compute generator matrix
############################################################

option (LIBLDPC_BUILD_GENERATOR "Build the ldpc_compute_generator binary (requires NTL)" OFF)

if (LIBLDPC_BUILD_GENERATOR)
	message(STATUS "Optional ldpc_compute_generator binary will be build. Disable this with -DLIBLDPC_BUILD_GENERATOR=Off")

	find_package(NTL 9.3.0 REQUIRED MODULE)
	find_package(GMP REQUIRED) # required by NTL, but not added to their public interface
	
	
	add_executable(ldpc_compute_generator ldpc_compute_generator.cpp)
	target_link_libraries(ldpc_compute_generator NTL::NTL)
	install(TARGETS ldpc_compute_generator DESTINATION bin)
else()
	message(STATUS "Optional ldpc_compute_generator binary is not build. Enable this with -DLIBLDPC_BUILD_GENERATOR=On")
endif()
