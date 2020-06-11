cmake_minimum_required(VERSION 3.10)

project(CALIB)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")

find_library(BLITZ_LIB libblitz.a PATHS ./lib)
find_library(FFTW3_LIB libfftw3.a)

if (NOT BLITZ_LIB)
message(FATAL_ERROR "File 'libblitz.a' not found! Blitz library must be installed!")
endif()
if (NOT FFTW3_LIB)
message(FATAL_ERROR "File 'libfftw3.a' not found! FFTW3 library must be installed!")
endif()

file(GLOB CALIB_SRC "./src/calib/*.cpp")
file(GLOB DEMO_SRC "./src/demo/*.cpp")
file(GLOB TEST_SRC "./src/test/*.cpp")

#Add target library

add_library(calib ${CALIB_SRC})

target_include_directories(calib PRIVATE ./inc)

target_link_libraries(
calib
${BLITZ_LIB} ${FFTW3_LIB}
)


#Add demo examples

set(DEMO)

foreach(demo_src ${DEMO_SRC})

get_filename_component(demo_bin ${demo_src} NAME)
string(REPLACE ".cpp" "" demo_bin ${demo_bin})

set(DEMO ${DEMO} ${demo_bin})

add_executable(${demo_bin} ${demo_src})

set_target_properties(
${demo_bin}
PROPERTIES RUNTIME_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/bin/demo
)
target_include_directories(
${demo_bin} PRIVATE ./inc
)
target_link_libraries(
${demo_bin} calib
)

endforeach()

add_custom_target(
demo DEPENDS
${DEMO}
)


#Add unit-tests

set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS unit_test_framework)

if (Boost_FOUND)

set(TEST)

foreach(test_src ${TEST_SRC})

get_filename_component(test_bin ${test_src} NAME)
string(REPLACE ".cpp" "" test_bin ${test_bin})

set(TEST ${TEST} ${test_bin})

add_executable(${test_bin} ${test_src})

set_target_properties(
${test_bin}
PROPERTIES RUNTIME_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/bin/test
)
target_include_directories(
${test_bin} PRIVATE ./inc
)
target_link_libraries(
${test_bin}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
target_link_libraries(
${test_bin} calib
)

endforeach()

add_custom_target(
test DEPENDS
${TEST}
)

endif()
