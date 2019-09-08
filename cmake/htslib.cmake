set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/contrib/htslib-install)
set(htslib_SRCDIR ${CMAKE_BINARY_DIR}/external/htslib)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl -lz")

ExternalProject_Add(htslib
        PREFIX ${htslib_PREFIX}
        SOURCE_DIR ${htslib_SRCDIR}
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf && ${htslib_SRCDIR}/configure --prefix=${htslib_INSTALL}
        BUILD_COMMAND make lib-static
        INSTALL_COMMAND make install prefix=${htslib_INSTALL}
        )

include_directories(${htslib_INSTALL}/include)
set(htslib_LIB ${htslib_INSTALL}/lib/libhts.a)