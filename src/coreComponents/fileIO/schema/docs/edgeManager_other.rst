

======================= =========================================================================================================== ========================================================= 
Name                    Type                                                                                                        Description                                               
======================= =========================================================================================================== ========================================================= 
domainBoundaryIndicator integer_array                                                                                               (no description available)                                
faceList                InterObjectRelation< Array< SortedArray< long, long >, 1, long, ChaiVector< SortedArray< long, long > > > > (no description available)                                
ghostRank               integer_array                                                                                               (no description available)                                
globalToLocalMap        map< long long, long, less< long long >, allocator< pair< long long const, long > > >                       (no description available)                                
isExternal              integer_array                                                                                               (no description available)                                
localToGlobalMap        globalIndex_array                                                                                           Array that contains a map from localIndex to globalIndex. 
nodeList                InterObjectRelation< Array< long, 2, long, ChaiVector< long > > >                                           (no description available)                                
parentIndex             localIndex_array                                                                                            Parent index of the edge.                                 
neighborData            node                                                                                                        :ref:`DATASTRUCTURE_neighborData`                         
sets                    node                                                                                                        :ref:`DATASTRUCTURE_sets`                                 
======================= =========================================================================================================== ========================================================= 


