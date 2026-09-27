add_executable(unit_tests
  test_toolchain_smoke.cpp
  test_eos.cpp
  test_mesh.cpp
)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain nlohmann_json::nlohmann_json eos mesh)

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
catch_discover_tests(unit_tests)
