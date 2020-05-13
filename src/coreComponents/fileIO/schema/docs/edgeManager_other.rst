

================================ ================================================================================================================ ==================================================================== 
Name                             Type                                                                                                             Description                                                          
================================ ================================================================================================================ ==================================================================== 
SIF_I                            real64_array                                                                                                     SIF_I of the edge.                                                   
SIF_II                           real64_array                                                                                                     SIF_II of the edge.                                                  
SIF_III                          real64_array                                                                                                     SIF_III of the edge.                                                 
childIndex                       localIndex_array                                                                                                 Child index of the edge.                                             
domainBoundaryIndicator          integer_array                                                                                                    (no description available)                                           
edgesToFractureConnectors        geosx_mapBase< long, long, std_integral_constant< bool, true > >                                                 A map of edge local indices to the fracture connector local indices. 
faceList                         geosx_InterObjectRelation< LvArray_ArrayOfSets< long, long > >                                                   (no description available)                                           
fractureConnectorsToEdges        localIndex_array                                                                                                 A map of fracture connector local indices to edge local indices.     
fractureConnectorsToElementIndex LvArray_ArrayOfArrays< long, long >                                                                              A map of fracture connector local indices face element local indices 
ghostRank                        integer_array                                                                                                    (no description available)                                           
globalToLocalMap                 geosx_mapBase< long long, long, std_integral_constant< bool, false > >                                           (no description available)                                           
isExternal                       integer_array                                                                                                    (no description available)                                           
localToGlobalMap                 globalIndex_array                                                                                                Array that contains a map from localIndex to globalIndex.            
nodeList                         geosx_InterObjectRelation< LvArray_Array< long, 2, camp_int_seq< long, 0l, 1l >, long, LvArray_NewChaiBuffer > > (no description available)                                           
parentIndex                      localIndex_array                                                                                                 Parent index of the edge.                                            
neighborData                     node                                                                                                             :ref:`DATASTRUCTURE_neighborData`                                    
sets                             node                                                                                                             :ref:`DATASTRUCTURE_sets`                                            
================================ ================================================================================================================ ==================================================================== 


