

================================ =========================================================================================================== ==================================================================== 
Name                             Type                                                                                                        Description                                                          
================================ =========================================================================================================== ==================================================================== 
domainBoundaryIndicator          integer_array                                                                                               (no description available)                                           
edgesToFractureConnectors        mapBase< long, long, integral_constant< bool, true > >                                                      A map of edge local indices to the fracture connector local indices. 
faceList                         InterObjectRelation< Array< SortedArray< long, long >, 1, long, ChaiVector< SortedArray< long, long > > > > (no description available)                                           
fractureConnectorsToEdges        localIndex_array                                                                                            A map of fracture connector local indices to edge local indices.     
fractureConnectorsToElementIndex ArrayOfArrays< long, long >                                                                                 A map of fracture connector local indices face element local indices 
ghostRank                        integer_array                                                                                               (no description available)                                           
globalToLocalMap                 mapBase< long long, long, integral_constant< bool, false > >                                                (no description available)                                           
isExternal                       integer_array                                                                                               (no description available)                                           
localToGlobalMap                 globalIndex_array                                                                                           Array that contains a map from localIndex to globalIndex.            
nodeList                         InterObjectRelation< Array< long, 2, long, ChaiVector< long > > >                                           (no description available)                                           
parentIndex                      localIndex_array                                                                                            Parent index of the edge.                                            
neighborData                     node                                                                                                        :ref:`DATASTRUCTURE_neighborData`                                    
sets                             node                                                                                                        :ref:`DATASTRUCTURE_sets`                                            
================================ =========================================================================================================== ==================================================================== 


