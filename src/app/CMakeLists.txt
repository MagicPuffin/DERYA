add_library(app STATIC driver.cpp)
target_include_directories(app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(app PUBLIC io hydro PRIVATE mesh eos timestep)

add_executable(hydro_run main.cpp)
target_link_libraries(hydro_run PRIVATE app io eos)
