find_package(PkgConfig)

pkg_check_modules(GMP QUIET GMP)

find_path(GMP_INCLUDE_DIR
    NAMES gmp.h
    PATHS ${PC_GMP_INCLUDE_DIRS}
    PATH_SUFFIXES GMP
)

find_library(GMP_LIBRARY
  NAMES gmp
  PATHS ${PC_GMP_LIBRARY_DIRS}
)


mark_as_advanced(GMP_FOUND GMP_INCLUDE_DIR GMP_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    FOUND_VAR GMP_FOUND
    REQUIRED_VARS
        GMP_LIBRARY
        GMP_INCLUDE_DIR
)

if(GMP_FOUND)
    set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    set(GMP_LIBRARIES ${GMP_LIBRARY})

    if(NOT TARGET GMP::GMP)
        add_library(GMP::GMP UNKNOWN IMPORTED)
    endif()

    set_target_properties(GMP::GMP PROPERTIES
        IMPORTED_LOCATION "${GMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
    )
endif()
