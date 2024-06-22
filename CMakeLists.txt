cmake_minimum_required(VERSION 3.5.0)
project(factoring_algorithms VERSION 0.1.0 LANGUAGES C CXX)

# set C++ flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
link_directories( $ENV{LD_LIBRARY_PATH} )


# setup ExternalProject and create an install location for it
include(ExternalProject)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

# Block Lanczos
ExternalProject_Add(extern-blanczos
    GIT_REPOSITORY https://github.com/SebWouters/blanczos
    PATCH_COMMAND bash -c "if ! grep -q PATCHED CMakeLists.txt$<SEMICOLON> then printf \"# PATCHED\\ninstall(TARGETS blanczos-shared)\\ninstall(FILES src/blanczos.h DESTINATION include)\\n\" >> CMakeLists.txt$<SEMICOLON> sed -i \"s\/std::cout\/$<SEMICOLON>\\\\/\\\\/std::cout/g\" src/blanczos.cpp$<SEMICOLON> fi"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

# Polynomial based on Eigen3 - non utile
# ExternalProject_Add(extern-polynomial
#     GIT_REPOSITORY https://github.com/jonathanventura/polynomial
#     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
# )

# googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  DOWNLOAD_EXTRACT_TIMESTAMP true
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# obtain lib and test soures
file(GLOB_RECURSE LIB_SRCS "src/*.cpp")
file(GLOB_RECURSE TEST_SRCS "test/*.cpp")
#include_directories(/usr/include/singular/)

# base library compiled from src/
add_library(factoring_algorithms_lib ${LIB_SRCS})
add_dependencies(factoring_algorithms_lib extern-blanczos)
target_link_libraries(factoring_algorithms_lib gmp)
target_link_libraries(factoring_algorithms_lib gmpxx)
target_link_libraries(factoring_algorithms_lib mpfr)
#target_link_libraries(factoring_algorithms_lib Singular)
target_link_libraries(factoring_algorithms_lib blanczos)
#target_link_libraries(factoring_algorithms_lib fflas)
#target_link_libraries(factoring_algorithms_lib givaro)
#target_link_libraries(factoring_algorithms_lib linbox)
target_link_libraries(factoring_algorithms_lib ntl)


# executable factoring_algorithms compiled from
# main.cpp and linked with the base library
add_executable(factoring_algorithms main.cpp)
target_link_libraries(factoring_algorithms factoring_algorithms_lib)


# executable all_tests to run tests compiled from all
# files in test/ and linked with the base library
enable_testing()
include(GoogleTest)
add_executable(all_tests ${TEST_SRCS})
target_link_libraries(all_tests factoring_algorithms_lib)
target_link_libraries(all_tests GTest::gtest_main)
gtest_discover_tests(all_tests)
