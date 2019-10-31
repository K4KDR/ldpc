find_package(PkgConfig)

pkg_check_modules(NTL QUIET NTL)

find_path(NTL_INCLUDE_DIR
    NAMES GF2.h
    PATHS ${PC_NTL_INCLUDE_DIRS}
    PATH_SUFFIXES NTL
)

find_library(NTL_LIBRARY
  NAMES ntl
  PATHS ${PC_NTL_LIBRARY_DIRS}
)


if ("${PC_NTL_VERSION}" STREQUAL "")
    # extract line with '#define NTL_VERSION "a.b.c"' from NTL/version.h
    file(STRINGS
         "${NTL_INCLUDE_DIR}/version.h"
         NTL_HEADER_VERSION_LINE
         REGEX "^[ \t]*#define[ \t]+NTL_VERSION[ \t]+\"[0-9.]+\"[ \t]*$"
    )

    # extract pure version number from line
    string(REGEX
           REPLACE
               "^[ \t]*#define[ \t]+NTL_VERSION[ \t]+\"([0-9.]+)\"[ \t]*$"
               "\\1"
           NTL_HEADER_VERSION ${NTL_HEADER_VERSION_LINE}
    )
    #message(STATUS "Found line '${NTL_HEADER_VERSION_LINE}' in version.h with version '${NTL_HEADER_VERSION}'")
    set(NTL_VERSION "${NTL_HEADER_VERSION}")
    if ("${NTL_VERSION}" STREQUAL "")
        message(ERROR "Cannot extract NTL version from NTL/version.h")
    else ()
        #message(STATUS "Found NTL header with version ${NTL_VERSION}")
        set(NTL_VERSION ${NTL_VERSION})
    endif ()
endif ()

set(NTL_VERSION ${NTL_VERSION})

mark_as_advanced(NTL_FOUND NTL_INCLUDE_DIR NTL_VERSION NTL_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NTL
    FOUND_VAR NTL_FOUND
    REQUIRED_VARS
        NTL_LIBRARY
        NTL_INCLUDE_DIR
    VERSION_VAR NTL_VERSION
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

find_package(GMP REQUIRED)

if(NTL_FOUND)
    set(NTL_INCLUDE_DIRS ${NTL_INCLUDE_DIR})
    set(NTL_LIBRARIES ${NTL_LIBRARY})

    if(NOT TARGET NTL::NTL)
        add_library(NTL::NTL UNKNOWN IMPORTED)
    endif()

    set_target_properties(NTL::NTL PROPERTIES
        IMPORTED_LOCATION "${NTL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NTL_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "Threads::Threads;GMP::GMP"
    )
endif()
