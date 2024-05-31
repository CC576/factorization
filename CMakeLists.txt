cmake_minimum_required(VERSION 3.5.0)
project(factoring_algorithms VERSION 0.1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 23)

include(ExternalProject)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(extern-blanczos
    GIT_REPOSITORY https://github.com/SebWouters/blanczos
    PATCH_COMMAND bash -c "if ! grep -q PATCHED CMakeLists.txt$<SEMICOLON> then printf \"# PATCHED\\ninstall(TARGETS blanczos-shared)\\ninstall(FILES src/blanczos.h DESTINATION include)\\n\" >> CMakeLists.txt$<SEMICOLON> fi"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

# non utile
# ExternalProject_Add(extern-polynomial
#     GIT_REPOSITORY https://github.com/jonathanventura/polynomial
#     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
# )

include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
#include_directories(/usr/include/singular/)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_executable(factoring_algorithms main.cpp)
target_link_libraries(factoring_algorithms gmp)
target_link_libraries(factoring_algorithms gmpxx)
target_link_libraries(factoring_algorithms mpfr)
#target_link_libraries(factoring_algorithms Singular)
target_link_libraries(factoring_algorithms blanczos)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
