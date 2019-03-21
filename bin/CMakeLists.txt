cmake_minimum_required(VERSION 3.0)

project(ldpc-bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/bin/cmake/Modules/")

############################################################
# Compute generator matrix
############################################################

option (BUILD_GENERATOR "Build the ldpc_compute_generator binary (requires NTL)" OFF)

if (BUILD_GENERATOR)
	message(STATUS "Optional ldpc_compute_generator binary will be build. Disable this with -DBUILD_GENERATOR=Off")

	find_package(NTL 9.3.0 REQUIRED MODULE)

	add_executable(ldpc_compute_generator ldpc_compute_generator.cpp)
	target_link_libraries(ldpc_compute_generator NTL::NTL)
	install(TARGETS ldpc_compute_generator DESTINATION bin)
else()
	message(STATUS "Optional ldpc_compute_generator binary is not build. Enable this with -DBUILD_GENERATOR=On")
endif()
