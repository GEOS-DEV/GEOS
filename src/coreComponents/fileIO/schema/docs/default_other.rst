

================================ ======================================================================================================================================= ============================================================= 
Name                             Type                                                                                                                                    Description                                                   
================================ ======================================================================================================================================= ============================================================= 
domainBoundaryIndicator          integer_array                                                                                                                           (no description available)                                    
edgeList                         InterObjectRelation< Array< Array< long, 1, int_seq< long, 0l >, long, NewChaiBuffer >, 1, int_seq< long, 0l >, long, NewChaiBuffer > > Map to the edges attached to each FaceElement.                
elementAperture                  real64_array                                                                                                                            The aperture of each FaceElement.                             
elementArea                      real64_array                                                                                                                            The area of each FaceElement.                                 
elementCenter                    r1_array                                                                                                                                The center of each FaceElement.                               
elementVolume                    real64_array                                                                                                                            The volume of each FaceElement.                               
faceList                         InterObjectRelation< Array< long, 2, int_seq< long, 0l, 1l >, long, NewChaiBuffer > >                                                   Map to the faces attached to each FaceElement.                
fractureElementsToCellIndices    localIndex_array2d                                                                                                                      A map of face element local indices to the cell local indices 
fractureElementsToCellRegions    localIndex_array2d                                                                                                                      A map of face element local indices to the cell local indices 
fractureElementsToCellSubRegions localIndex_array2d                                                                                                                      A map of face element local indices to the cell local indices 
ghostRank                        integer_array                                                                                                                           (no description available)                                    
globalToLocalMap                 mapBase< long long, long, integral_constant< bool, false > >                                                                            (no description available)                                    
isExternal                       integer_array                                                                                                                           (no description available)                                    
localToGlobalMap                 globalIndex_array                                                                                                                       Array that contains a map from localIndex to globalIndex.     
nodeList                         InterObjectRelation< Array< Array< long, 1, int_seq< long, 0l >, long, NewChaiBuffer >, 1, int_seq< long, 0l >, long, NewChaiBuffer > > Map to the nodes attached to each FaceElement.                
ConstitutiveModels               node                                                                                                                                    :ref:`DATASTRUCTURE_ConstitutiveModels`                       
neighborData                     node                                                                                                                                    :ref:`DATASTRUCTURE_neighborData`                             
sets                             node                                                                                                                                    :ref:`DATASTRUCTURE_sets`                                     
================================ ======================================================================================================================================= ============================================================= 


