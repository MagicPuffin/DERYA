add_executable(unit_tests
  test_toolchain_smoke.cpp
  test_eos.cpp
)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain nlohmann_json::nlohmann_json eos)

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
catch_discover_tests(unit_tests)
