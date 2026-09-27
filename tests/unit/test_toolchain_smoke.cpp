// Task 0 smoke test: confirms the build toolchain itself works before any
// module depends on it. Two things checked:
//   1. Catch2 is correctly fetched/linked (this test running at all proves it)
//   2. nlohmann::json is correctly fetched/linked (deck format dependency)
// This file should stay this trivial permanently — it is not the place for
// real tests once modules exist.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

TEST_CASE("toolchain: Catch2 runs", "[smoke]") {
    REQUIRE(1 + 1 == 2);
}

TEST_CASE("toolchain: nlohmann::json is linked and functional", "[smoke]") {
    nlohmann::json j;
    j["end_time"] = 0.2;
    j["mesh"]["nx"] = 100;
    REQUIRE(j["end_time"].get<double>() == 0.2);
    REQUIRE(j["mesh"]["nx"].get<int>() == 100);
    REQUIRE(j.dump() == R"({"end_time":0.2,"mesh":{"nx":100}})");
}
