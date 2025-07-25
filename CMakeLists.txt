cmake_minimum_required(VERSION 3.18)
project(impaqt)

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files
file(GLOB IMPAQT_FILES
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

# Include paths
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/ext/bamtools/src
    ${PROJECT_SOURCE_DIR}/ext/seqan/include
    ${GTEST_INCLUDE_DIRS}
)

# BamTools: Prevent installation
set(BAMTOOLS_SUBPROJECT ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/bamtools/src/api EXCLUDE_FROM_ALL)

# Executable
add_executable(impaqt ${IMPAQT_FILES})
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
