cmake_minimum_required(VERSION "3.18.0")
project(impaqt)

file(GLOB IMPACT_FILES
    ${PROJECT_SOURCE_DIR}/src/impaqt.cpp
)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -pthread -std=c++14")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -pthread -std=c++14 -ggdb3")

# Dependencies
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

include(cmake/gtest.cmake)
message( status "include dirs: " ${GTEST_INCLUDE_DIRS})

# Statically link program and depedencies
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/ext/bamtools/src")
include_directories("${PROJECT_SOURCE_DIR}/ext/seqan/include")
include_directories("${GTEST_INCLUDE_DIRS} ${GTEST_INCLUDE_DIR}")

add_definitions(${SEQAN_DEFINITIONS})
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/bamtools/src/api)

enable_testing()

# Program
add_executable(impaqt ${IMPACT_FILES})
target_link_libraries (impaqt ${SEQAN_LIBRARIES})
target_link_libraries(impaqt ${ZLIB_LIBRARIES})
target_link_libraries(impaqt BamTools)

# Test
link_directories(${GTEST_LIBS_DIR})
add_executable(impaqt_test test/impaqt_test.cpp)
target_link_libraries(impaqt_test gtest gtest_main)
target_link_libraries(impaqt_test BamTools)
add_test(impaqt_test impaqt_test)


# Install the executable
install(TARGETS impaqt DESTINATION bin)
