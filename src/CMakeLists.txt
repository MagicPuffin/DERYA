# Each module below is its own CMake target with an explicit public include
# directory. Modules may only link against the targets listed in their own
# target_link_libraries() — a module reaching into another module's private
# headers is a build failure, not a convention. See docs/DESIGN.md.
#
# Populated incrementally per the v0.1 backlog (Task 1 onward). Empty
# add_subdirectory calls are commented out until each module has real content,
# so the tree structure is visible from day one without breaking the build.

add_subdirectory(eos)
add_subdirectory(mesh)
# add_subdirectory(bc)
# add_subdirectory(timestep)
# add_subdirectory(hydro)
# add_subdirectory(sync)
# add_subdirectory(fields)
# add_subdirectory(io)
# add_subdirectory(app)
