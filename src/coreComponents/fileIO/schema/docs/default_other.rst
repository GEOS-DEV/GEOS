

================================ ================================================================================================================ ========================================================= 
Name                             Type                                                                                                             Description                                               
================================ ================================================================================================================ ========================================================= 
domainBoundaryIndicator          integer_array                                                                                                    (no description available)                                
edgeList                         geosx_InterObjectRelation< LvArray_ArrayOfArrays< long, long > >                                                 Map to the edges.                                         
elementAperture                  real64_array                                                                                                     The aperture of each EmbeddedSurface.                     
elementArea                      real64_array                                                                                                     The area of each EmbeddedSurface element.                 
elementCenter                    r1_array                                                                                                         The center of each EmbeddedSurface element.               
elementVolume                    real64_array                                                                                                     The volume of each EmbeddedSurface element.               
fractureElementsToCellIndices    localIndex_array                                                                                                 Map to the cells.                                         
fractureElementsToRegionIndex    localIndex_array                                                                                                 Map to the region cut by each EmbeddedSurface.            
fractureElementsToSubRegionIndex localIndex_array                                                                                                 Map to the subregion cut by each EmbeddedSurface.         
ghostRank                        integer_array                                                                                                    (no description available)                                
globalToLocalMap                 geosx_mapBase< long long, long, std_integral_constant< bool, false > >                                           (no description available)                                
isExternal                       integer_array                                                                                                    (no description available)                                
localToGlobalMap                 globalIndex_array                                                                                                Array that contains a map from localIndex to globalIndex. 
nodeList                         geosx_InterObjectRelation< LvArray_Array< long, 2, camp_int_seq< long, 0l, 1l >, long, LvArray_NewChaiBuffer > > Map to the nodes attached to each EmbeddedSurface.        
normalVector                     r1_array                                                                                                         Unit normal vector to the embedded surface.               
ConstitutiveModels               node                                                                                                             :ref:`DATASTRUCTURE_ConstitutiveModels`                   
neighborData                     node                                                                                                             :ref:`DATASTRUCTURE_neighborData`                         
sets                             node                                                                                                             :ref:`DATASTRUCTURE_sets`                                 
================================ ================================================================================================================ ========================================================= 


