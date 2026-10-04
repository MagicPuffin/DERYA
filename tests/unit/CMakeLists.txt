add_executable(unit_tests
  test_toolchain_smoke.cpp
  test_eos.cpp
  test_mesh.cpp
  test_bc.cpp
  test_hydro.cpp
  test_timestep.cpp
)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain nlohmann_json::nlohmann_json eos mesh bc hydro sync timestep)
target_compile_definitions(unit_tests PRIVATE
  ORACLE_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/oracle_data")

include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
catch_discover_tests(unit_tests)
