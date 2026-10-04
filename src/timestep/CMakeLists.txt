add_library(timestep STATIC timestep.cpp)
target_include_directories(timestep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(timestep PUBLIC hydro eos mesh PRIVATE sync)
