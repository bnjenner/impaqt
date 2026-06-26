cmake_minimum_required(VERSION 3.18)
project(impaqt)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Warning flags for first-party targets only (vendored ext/ deps are included
# as SYSTEM headers and built separately, so they are not affected).
set(IMPAQT_WARNINGS -Wall -Wextra)

# Source files (listed explicitly; GLOB does not re-run when files are added)
set(IMPAQT_FILES
	${PROJECT_SOURCE_DIR}/src/impaqt.cpp
	${PROJECT_SOURCE_DIR}/src/AnnotationList.cpp
	${PROJECT_SOURCE_DIR}/src/ClusterList.cpp
	${PROJECT_SOURCE_DIR}/src/DBSCAN.cpp
	${PROJECT_SOURCE_DIR}/src/ContainmentList.cpp
	${PROJECT_SOURCE_DIR}/src/AssignClusters.cpp
	${PROJECT_SOURCE_DIR}/src/utils.cpp
)

# Dependencies
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Build options
#  USE_SYSTEM_BAMTOOLS  link an externally-provided bamtools (found via pkg-config,
#                       e.g. from conda or the system) instead of fetching+building
#                       it. Use for packaging (conda/bioconda) where deps come from
#                       the package manager and the build must not download anything.
#  IMPAQT_BUILD_TESTS   build the unit tests (fetches GoogleTest). Turn OFF for
#                       packaging so no test deps are needed at build time.
option(USE_SYSTEM_BAMTOOLS "Link an external bamtools instead of fetching it" OFF)
option(IMPAQT_BUILD_TESTS  "Build the unit tests (fetches GoogleTest)"        ON)

include(FetchContent)

# --- BamTools ---
if(USE_SYSTEM_BAMTOOLS)
    # Find an installed bamtools via its pkg-config file (bamtools-1.pc) and expose
    # it under the same 'BamTools' target name the rest of this file links against.
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(BAMTOOLS REQUIRED IMPORTED_TARGET bamtools-1)
    add_library(BamTools INTERFACE)
    # A system bamtools is usually a static lib whose zlib dependency is only a
    # pkg-config *private* requirement (not pulled in for non-static linking), so
    # link zlib explicitly or its BGZF deflate/inflate/crc32 symbols stay undefined.
    target_link_libraries(BamTools INTERFACE PkgConfig::BAMTOOLS ZLIB::ZLIB)
    set(BAMTOOLS_INCLUDE_DIR ${BAMTOOLS_INCLUDE_DIRS})
else()
    # Fetch bamtools v2.5.3, pinned to an immutable commit SHA so builds are
    # reproducible. v2.5.3 already carries the std::-qualified C-function fix the
    # old vendored copy was hand-patched for, so it builds clean under C++17. The
    # BamTools library target lives in src/ in 2.5.3, so use a full MakeAvailable.
    FetchContent_Declare(
        bamtools
        GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
        GIT_TAG 7e5c183cdf37a6ab216402b26c07296807416472  # v2.5.3
        SYSTEM            # treat bamtools headers as system (suppress its warnings)
        EXCLUDE_FROM_ALL  # build only the BamTools lib we link; skip toolkit & install
    )
    FetchContent_MakeAvailable(bamtools)
    # bamtools registers 2 ctest cases (bamtools_help/_stats) that invoke its toolkit
    # binary, which we don't build (EXCLUDE_FROM_ALL); they're in a subdirectory scope
    # (unreachable by set_tests_properties), so tell ctest to ignore them.
    file(WRITE "${CMAKE_BINARY_DIR}/CTestCustom.cmake"
         "set(CTEST_CUSTOM_TESTS_IGNORE bamtools_help bamtools_stats)\n")
    set(BAMTOOLS_INCLUDE_DIR ${bamtools_SOURCE_DIR}/src)
endif()

# --- GoogleTest (unit tests only) ---
if(IMPAQT_BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG 58d77fa8070e8cec2dc1ed015d66b454c8d78850  # release-1.12.1
        SYSTEM
        EXCLUDE_FROM_ALL  # don't add gtest's targets/install rules to our build & install
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)  # don't override the parent's CRT (Windows)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)           # don't install gtest with impaqt
    FetchContent_MakeAvailable(googletest)
endif()

# Include paths
include_directories(SYSTEM
    ${PROJECT_SOURCE_DIR}/include
    ${BAMTOOLS_INCLUDE_DIR}
)

# Executable
add_executable(impaqt ${IMPAQT_FILES})
target_compile_options(impaqt PRIVATE ${IMPAQT_WARNINGS})
target_link_libraries(impaqt
    BamTools
    ${ZLIB_LIBRARIES}
    Threads::Threads
)

# Install the main binary only
install(TARGETS impaqt DESTINATION bin)

# -----------------------
# Test Section
# -----------------------

if(IMPAQT_BUILD_TESTS)

enable_testing()


# Test: annotation_test
add_executable(annotation_test
    ${PROJECT_SOURCE_DIR}/test/AnnotationList_test.cpp
)
target_sources(annotation_test
    PRIVATE ${PROJECT_SOURCE_DIR}/src/AnnotationList.cpp
    ${PROJECT_SOURCE_DIR}/src/ClusterList.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
)
target_link_libraries(annotation_test
    gtest gtest_main
    BamTools
)
add_test(NAME annotation_test COMMAND annotation_test)


# Test: cluster_test
add_executable(cluster_test
    ${PROJECT_SOURCE_DIR}/test/ClusterList_test.cpp
)
target_sources(cluster_test
    PRIVATE ${PROJECT_SOURCE_DIR}/src/ClusterList.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
)
target_link_libraries(cluster_test
    gtest gtest_main
    BamTools
)
add_test(NAME cluster_test COMMAND cluster_test)


# Test: dbscan_test
add_executable(dbscan_test
    ${PROJECT_SOURCE_DIR}/test/DBSCAN_test.cpp
)
target_sources(dbscan_test
    PRIVATE ${PROJECT_SOURCE_DIR}/src/ClusterList.cpp
    ${PROJECT_SOURCE_DIR}/src/ContainmentList.cpp
    ${PROJECT_SOURCE_DIR}/src/DBSCAN.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
)
target_link_libraries(dbscan_test
    gtest gtest_main
    BamTools
)
add_test(NAME dbscan_test COMMAND dbscan_test)


# Test: assign_test
add_executable(assign_test
    ${PROJECT_SOURCE_DIR}/test/AssignClusters_test.cpp
)
target_sources(assign_test
    PRIVATE ${PROJECT_SOURCE_DIR}/src/AssignClusters.cpp
    ${PROJECT_SOURCE_DIR}/src/ClusterList.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
)
target_link_libraries(assign_test
    gtest gtest_main
    BamTools
)
add_test(NAME assign_test COMMAND assign_test)

endif()  # IMPAQT_BUILD_TESTS
