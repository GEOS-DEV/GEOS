

================================ ================================================================================== ========================================================= 
Name                             Type                                                                               Description                                               
================================ ================================================================================== ========================================================= 
domainBoundaryIndicator          Array< int, 1, int_seq< long, 0l >, long, ChaiBuffer >                             (no description available)                                
edgeList                         InterObjectRelation< ArrayOfArrays< long, long, ChaiBuffer > >                     Map to the edges.                                         
elementAperture                  Array< double, 1, int_seq< long, 0l >, long, ChaiBuffer >                          The aperture of each EmbeddedSurface.                     
elementArea                      Array< double, 1, int_seq< long, 0l >, long, ChaiBuffer >                          The area of each EmbeddedSurface element.                 
elementCenter                    Array< double, 2, int_seq< long, 0l, 1l >, long, ChaiBuffer >                      (no description available)                                
elementVolume                    Array< double, 1, int_seq< long, 0l >, long, ChaiBuffer >                          (no description available)                                
fractureElementsToCellIndices    Array< long, 1, int_seq< long, 0l >, long, ChaiBuffer >                            Map to the cells.                                         
fractureElementsToRegionIndex    Array< long, 1, int_seq< long, 0l >, long, ChaiBuffer >                            Map to the region cut by each EmbeddedSurface.            
fractureElementsToSubRegionIndex Array< long, 1, int_seq< long, 0l >, long, ChaiBuffer >                            Map to the subregion cut by each EmbeddedSurface.         
ghostRank                        Array< int, 1, int_seq< long, 0l >, long, ChaiBuffer >                             (no description available)                                
globalToLocalMap                 mapBase< long long, long, integral_constant< bool, false > >                       (no description available)                                
isExternal                       Array< int, 1, int_seq< long, 0l >, long, ChaiBuffer >                             (no description available)                                
localToGlobalMap                 Array< long long, 1, int_seq< long, 0l >, long, ChaiBuffer >                       Array that contains a map from localIndex to globalIndex. 
nodeList                         InterObjectRelation< Array< long, 2, int_seq< long, 0l, 1l >, long, ChaiBuffer > > Map to the nodes attached to each EmbeddedSurface.        
normalVector                     Array< R1TensorT<3>, 1, int_seq< long, 0l >, long, ChaiBuffer >                    Unit normal vector to the embedded surface.               
numEdgesPerElement               localIndex                                                                         (no description available)                                
numFacesPerElement               localIndex                                                                         (no description available)                                
numNodesPerElement               localIndex                                                                         (no description available)                                
ConstitutiveModels               node                                                                               :ref:`DATASTRUCTURE_ConstitutiveModels`                   
neighborData                     node                                                                               :ref:`DATASTRUCTURE_neighborData`                         
sets                             node                                                                               :ref:`DATASTRUCTURE_sets`                                 
================================ ================================================================================== ========================================================= 


