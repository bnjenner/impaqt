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
    ${PROJECT_SOURCE_DIR}/lib
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

# Test: impaqt_test
#add_executable(impaqt_test test/impaqt_test.cpp)
#target_link_libraries(impaqt_test
#    gtest gtest_main
#    BamTools
#)
#add_test(NAME impaqt_test COMMAND impaqt_test)

# Test: ClusterList_test
#add_executable(ClusterList_test test/ClusterList_test.cpp)
#target_link_libraries(ClusterList_test
#    gtest gtest_main
#    BamTools
#)
#add_test(NAME ClusterList_test COMMAND ClusterList_test)
