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

# GoogleTest
include(cmake/gtest.cmake)
message(STATUS "GTest include dirs: ${GTEST_INCLUDE_DIRS}")

# BamTools via FetchContent, pinned to v2.5.3 by immutable commit SHA.
# A full SHA (not a bare tag) can't be force-moved upstream, so builds are
# reproducible. v2.5.3 already carries the std::-qualified C-function fix that the
# previously-vendored copy was hand-patched for, so it builds clean under C++17.
include(FetchContent)
FetchContent_Declare(
    bamtools
    GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
    GIT_TAG 7e5c183cdf37a6ab216402b26c07296807416472  # v2.5.3
    SYSTEM  # treat bamtools headers as system (suppress its warnings under -Wall -Wextra)
)
# Full MakeAvailable: in 2.5.3 the BamTools library target lives in src/ (not
# src/api), so we build the project as upstream intends. It also builds the
# bamtools toolkit + registers 2 of its own ctest cases (bamtools_help/_stats);
# both pass and serve as a sanity check that the pinned dep itself works.
FetchContent_MakeAvailable(bamtools)

# Include paths
include_directories(SYSTEM
    ${PROJECT_SOURCE_DIR}/include
    ${bamtools_SOURCE_DIR}/src
    ${GTEST_INCLUDE_DIRS}
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
link_directories(${GTEST_LIBS_DIR})


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
