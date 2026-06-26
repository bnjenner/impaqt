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

# Third-party dependencies via FetchContent, each pinned to an immutable commit
# SHA (not a bare tag, which can be force-moved upstream) so builds are reproducible.
include(FetchContent)

# GoogleTest (release-1.12.1) -- tests only.
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

# BamTools (v2.5.3). v2.5.3 already carries the std::-qualified C-function fix the
# previously-vendored copy was hand-patched for, so it builds clean under C++17.
# In 2.5.3 the BamTools library target lives in src/ (not src/api), so we build the
# project as upstream intends via a full MakeAvailable; this also builds the toolkit
# and registers 2 of its own ctest cases (bamtools_help/_stats), which pass and
# sanity-check the pinned dep.
FetchContent_Declare(
    bamtools
    GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
    GIT_TAG 7e5c183cdf37a6ab216402b26c07296807416472  # v2.5.3
    SYSTEM            # treat bamtools headers as system (suppress its warnings under -Wall -Wextra)
    EXCLUDE_FROM_ALL  # build only the BamTools lib we link; skip its toolkit & install rules
)
FetchContent_MakeAvailable(bamtools)

# bamtools' root CMakeLists registers two ctest cases (bamtools_help/_stats) that
# invoke its toolkit binary, which we don't build (EXCLUDE_FROM_ALL). They are
# defined in a subdirectory scope (so set_tests_properties can't reach them from
# here); instead tell ctest to ignore them via CTestCustom, so `ctest` runs only
# impaqt's own suite.
file(WRITE "${CMAKE_BINARY_DIR}/CTestCustom.cmake"
     "set(CTEST_CUSTOM_TESTS_IGNORE bamtools_help bamtools_stats)\n")

# Include paths
include_directories(SYSTEM
    ${PROJECT_SOURCE_DIR}/include
    ${bamtools_SOURCE_DIR}/src
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
