set(hisat2_PREFIX ${CMAKE_BINARY_DIR}/contrib/hisat2-prefix)
set(hisat2_INSTALL ${CMAKE_BINARY_DIR}/contrib/hisat2-install)
set(hisat2_SRCDIR ${CMAKE_BINARY_DIR}/external/hisat2)

ExternalProject_Add(hisat2
        PREFIX ${hisat2_PREFIX}
        SOURCE_DIR ${hisat2_SRCDIR}
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND cmake -DCMAKE_BUILD_TYPE=Release ${hisat2_SRCDIR} && make -C ${hisat2_SRCDIR}
        INSTALL_COMMAND cmake -E echo "Skipping install step."
)