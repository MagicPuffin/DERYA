add_library(hydro STATIC hydro.cpp)
target_include_directories(hydro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hydro PUBLIC mesh bc eos)
