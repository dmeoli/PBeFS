# the library with the runtimes of its features; FastFlow, which no port
# provides and which would be downloaded at configure time, is left out
vcpkg_download_distfile(ARCHIVE
    URLS "https://github.com/dmeoli/PBeFS/releases/download/v${VERSION}/PBeFS-${VERSION}.tar.gz"
    FILENAME "PBeFS-${VERSION}.tar.gz"
    SHA512 0
)
vcpkg_extract_source_archive(SOURCE_PATH ARCHIVE "${ARCHIVE}")

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        openmp PBEFS_WITH_OPENMP
        mpi    PBEFS_WITH_MPI
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        ${FEATURE_OPTIONS}
        -DPBEFS_WITH_FASTFLOW=OFF
        -DPBEFS_FETCH_FASTFLOW=OFF
        -DPBEFS_BUILD_TESTS=OFF
        -DPBEFS_BUILD_EXAMPLES=OFF
        -DPBEFS_BUILD_TOOLS=OFF
)
vcpkg_cmake_install()

# the README and the LICENSE that the project installs in share/PBeFS,
# which on Windows is the directory share/pbefs of the port
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/share/PBeFS"
                    "${CURRENT_PACKAGES_DIR}/debug/share/PBeFS")
vcpkg_cmake_config_fixup(PACKAGE_NAME pbefs CONFIG_PATH lib/cmake/PBeFS)
vcpkg_fixup_pkgconfig()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include"
                    "${CURRENT_PACKAGES_DIR}/debug/share")
file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage"
     DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
