add_library(io STATIC deck.cpp dump.cpp)
target_include_directories(io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(io PUBLIC mesh bc eos hydro timestep PRIVATE nlohmann_json::nlohmann_json)
