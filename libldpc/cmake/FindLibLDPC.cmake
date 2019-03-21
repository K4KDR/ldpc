#FindLibLDPC.cmake
#
# Finds the LDPC library
#
# This will define the following variables
#
#    LIBLDPC_FOUND
#    LIBLDPC_INCLUDE_DIRS
#
# and the following imported targets
#
#     libldpc::libldpc

find_package(PkgConfig)
pkg_check_modules(PC_libldpc QUIET LibLDPC)

find_path(LIBLDPC_INCLUDE_DIR
    NAMES ldpc.h
    PATHS ${PC_LIBLDPC_INCLUDE_DIRS}
    PATH_SUFFIXES ldpc
)

set(LIBLDPC_VERSION ${PC_LIBLDPC_VERSION})

mark_as_advanced(LIBLDPC_FOUND LIBDPC_INCLUDE_DIR LIBLDPC_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBLDPC
    REQUIRED_VARS LIBLDPC_INCLUDE_DIR
    VERSION_VAR LIBLDPC_VERSION
)

if(LIBLDPC_FOUND)
    #Set include dirs to parent, to enable includes like #include <ldpc/encoder.h>
    get_filename_component(LIBLDPC_INCLUDE_DIRS ${LIBLDPC_INCLUDE_DIR} DIRECTORY)
endif()

if(LIBLDPC_FOUND AND NOT TARGET libldpc::libldpc)
    add_library(libldpc::libldpc INTERFACE IMPORTED)
    set_target_properties(libldpc::libldpc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBLDPC_INCLUDE_DIRS}"
    )
endif()
