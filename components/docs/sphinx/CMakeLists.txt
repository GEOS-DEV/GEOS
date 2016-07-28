
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

add_sphinx_target( sidre_docs )

if(DOXYGEN_FOUND)
   add_dependencies(sidre_docs sidre_doxygen) 
endif()

