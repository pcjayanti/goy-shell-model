cmake_minimum_required(VERSION 3.16)
project(goy_shell_model VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release by default: a Debug build of this integrator is ~20x slower, which
# makes any timing comparison against the 2016 C version meaningless.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

add_library(goy_warnings INTERFACE)
if(MSVC)
  target_compile_options(goy_warnings INTERFACE /W4)
else()
  # -Wconversion matters here: it is the warning that catches the integer
  # division in `pow(eta, -1/6)`, which -Wall -Wextra let through silently.
  target_compile_options(goy_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion)
endif()

add_library(goy_core
  src/model.cpp
  src/parameters.cpp
  src/recorder.cpp
)
target_include_directories(goy_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(goy_core PRIVATE goy_warnings)

find_package(Threads REQUIRED)

add_executable(goy src/main.cpp)
target_link_libraries(goy PRIVATE goy_core goy_warnings Threads::Threads)

enable_testing()
add_executable(goy_tests tests/test_goy.cpp)
target_link_libraries(goy_tests PRIVATE goy_core goy_warnings)
add_test(NAME goy_tests COMMAND goy_tests)

# Optional: -march=native gives a useful speedup but produces a binary that
# only runs on this machine. Off by default; enable with -DGOY_NATIVE=ON.
option(GOY_NATIVE "Optimise for the building machine's CPU" OFF)
if(GOY_NATIVE AND NOT MSVC)
  target_compile_options(goy_core PRIVATE -march=native)
  target_compile_options(goy PRIVATE -march=native)
endif()
