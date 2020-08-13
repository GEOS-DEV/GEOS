

================================ ==================================================================================== ==================================================================== 
Name                             Type                                                                                 Description                                                          
================================ ==================================================================================== ==================================================================== 
connectivityIndex                real64_array                                                                         Connectivity index of each EmbeddedSurface.                          
domainBoundaryIndicator          integer_array                                                                        (no description available)                                           
elementAperture                  real64_array                                                                         The aperture of each EmbeddedSurface.                                
elementArea                      real64_array                                                                         The area of each EmbeddedSurface element.                            
elementCenter                    real64_array2d                                                                       The center of each EmbeddedSurface element.                          
elementVolume                    real64_array                                                                         The volume of each EmbeddedSurface element.                          
fractureElementsToCellIndices    localIndex_array                                                                     Map to the cells.                                                    
fractureElementsToRegionIndex    localIndex_array                                                                     Map to the region cut by each EmbeddedSurface.                       
fractureElementsToSubRegionIndex localIndex_array                                                                     Map to the subregion cut by each EmbeddedSurface.                    
ghostRank                        integer_array                                                                        (no description available)                                           
globalToLocalMap                 geosx_mapBase< long long, long, std_integral_constant< bool, false > >               (no description available)                                           
isExternal                       integer_array                                                                        (no description available)                                           
localToGlobalMap                 globalIndex_array                                                                    Array that contains a map from localIndex to globalIndex.            
nodeList                         geosx_InterObjectRelation< LvArray_ArrayOfArrays< long, long, LvArray_ChaiBuffer > > Map to the nodes attached to each EmbeddedSurface.                   
normalVector                     r1_array                                                                             Unit normal vector to the embedded surface.                          
numEdgesPerElement               localIndex                                                                           (no description available)                                           
numFacesPerElement               localIndex                                                                           (no description available)                                           
numNodes                         integer_array                                                                        Number of nodes of each EmbeddedSurface.                             
numNodesPerElement               localIndex                                                                           (no description available)                                           
tangentVector1                   r1_array                                                                             Unit vector in the first tangent direction to the embedded surface.  
tangentVector2                   r1_array                                                                             Unit vector in the second tangent direction to the embedded surface. 
ConstitutiveModels               node                                                                                 :ref:`DATASTRUCTURE_ConstitutiveModels`                              
neighborData                     node                                                                                 :ref:`DATASTRUCTURE_neighborData`                                    
sets                             node                                                                                 :ref:`DATASTRUCTURE_sets`                                            
================================ ==================================================================================== ==================================================================== 


