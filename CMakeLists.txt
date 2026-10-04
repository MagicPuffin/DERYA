cmake_minimum_required(VERSION 3.20)
project(hydro LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_MPI "Build p4est/libsc against real MPI (OFF = single-rank stub build)" OFF)
option(BUILD_TESTS "Build the Catch2 unit test suite" ON)

# Optimized with debug symbols by default: the integration tests are ~12x
# slower unoptimized. Debug builds go in their own tree (the "debug" preset,
# build-debug/); see docs/decisions.md, "Build trees".
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

include(FetchContent)

# --- Catch2 (unit test framework, Tier 1) -----------------------------------
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1
)
FetchContent_MakeAvailable(Catch2)

# --- nlohmann::json (deck format) -------------------------------------------
FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
)
FetchContent_MakeAvailable(json)

# p4est/libsc intentionally NOT fetched here — deferred to v0.4 per the
# v0.1 backlog (v0.1-v0.3 use a minimal non-p4est structured mesh).

add_subdirectory(src)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/unit)
endif()
