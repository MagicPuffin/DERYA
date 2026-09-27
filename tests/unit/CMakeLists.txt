add_executable(unit_tests
  test_toolchain_smoke.cpp
)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain nlohmann_json::nlohmann_json)

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
catch_discover_tests(unit_tests)
