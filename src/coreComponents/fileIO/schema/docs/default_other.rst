

================================ ============================================================ ========================================================= 
Name                             Type                                                         Description                                               
================================ ============================================================ ========================================================= 
domainBoundaryIndicator          integer_array                                                (no description available)                                
elementAperture                  real64_array                                                 The aperture of each EmbeddedSurface.                     
elementArea                      real64_array                                                 The area of each EmbeddedSurface element.                 
elementCenter                    r1_array                                                     The center of each EmbeddedSurface element.               
elementVolume                    real64_array                                                 The volume of each EmbeddedSurface element.               
fractureElementsToCellIndices    localIndex_array                                             Map to the cells.                                         
fractureElementsToRegionIndex    localIndex_array                                             Map to the region cut by each EmbeddedSurface.            
fractureElementsToSubRegionIndex localIndex_array                                             Map to the subregion cut by each EmbeddedSurface.         
ghostRank                        integer_array                                                (no description available)                                
globalToLocalMap                 mapBase< long long, long, integral_constant< bool, false > > (no description available)                                
isExternal                       integer_array                                                (no description available)                                
localToGlobalMap                 globalIndex_array                                            Array that contains a map from localIndex to globalIndex. 
normalVector                     r1_array                                                     Unit normal vector to the embedded surface.               
ConstitutiveModels               node                                                         :ref:`DATASTRUCTURE_ConstitutiveModels`                   
neighborData                     node                                                         :ref:`DATASTRUCTURE_neighborData`                         
sets                             node                                                         :ref:`DATASTRUCTURE_sets`                                 
================================ ============================================================ ========================================================= 


