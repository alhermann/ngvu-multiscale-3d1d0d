add_executable(elementparametrization elementparametrization.cc)
target_link_libraries(elementparametrization ${DUNE_LIBS})

add_executable(parametrized-refinement parametrized-refinement.cc)
target_link_libraries(parametrized-refinement ${DUNE_LIBS})
