cmake_minimum_required(VERSION 3.0)

project(ldpc-lib VERSION 0.1.0 LANGUAGES CXX)

##############################################
# Create target and set properties

set(libldpc_SOURCES
	src/decoder.cpp
	src/encoder.cpp
	src/ldpc.cpp
)

add_library(libldpc SHARED ${libldpc_SOURCES})

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(libldpc::libldpc ALIAS libldpc)

target_include_directories(libldpc
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Set compiler to strict mode when compiling this binary
target_compile_options(libldpc
    PRIVATE
        $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
            -Werror -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion>
        $<$<C_COMPILER_ID:MSVC>:
            /WX /W4>
)

##############################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ldpc)

install(TARGETS libldpc
    EXPORT libldpc-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Export the targets to a script
install(EXPORT libldpc-targets
  FILE
    LibLDPCTargets.cmake
  NAMESPACE
    libldpc::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/LibLDPCConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/LibLDPCConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/LibLDPCConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(FILES
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindLibLDPC.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/LibLDPCConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/LibLDPCConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindLibLDPC.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FindLibLDPC.cmake
    COPYONLY)

export(EXPORT libldpc-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/LibLDPCTargets.cmake NAMESPACE libldpc::)

#Register package in user's package registry
export(PACKAGE libldpc)
