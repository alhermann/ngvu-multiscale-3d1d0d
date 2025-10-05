add_definitions(-DDUNE_FOAMGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

dune_add_test(SOURCES foamgrid-test.cc)
dune_add_test(SOURCES global-refine-test.cc)
dune_add_test(SOURCES local-refine-test.cc)
dune_add_test(SOURCES growth-test-1d.cc)
dune_add_test(SOURCES growth-test-2d.cc)
dune_add_test(SOURCES setposition-test.cc)
