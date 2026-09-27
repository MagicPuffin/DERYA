add_library(bc STATIC bc.cpp)
target_include_directories(bc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(bc PUBLIC mesh)
