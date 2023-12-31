cmake_minimum_required(VERSION 3.10)

project(ldpc LANGUAGES CXX VERSION 1.0)

set(CMAKE_MODULE_PATH $7CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

add_subdirectory(libldpc)
add_subdirectory(bin)

option(LDPC_UNITTESTS "Enable unittests" Off)

if(LDPC_UNITTESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Add uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

#Export the targets to a script
install(EXPORT ldpcTargets
  FILE ldpcTargets.cmake
  NAMESPACE ldpc::
  DESTINATION lib/cmake/ldpc
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ldpcConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

#Install the config, configversion and custom find modules
install(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ldpcConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/ldpcConfigVersion.cmake
    DESTINATION lib/cmake/ldpc
)

