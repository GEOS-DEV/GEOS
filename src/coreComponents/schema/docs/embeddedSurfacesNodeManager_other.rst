

======================= ================================================================================== ============================================= ================================================================ 
Name                    Type                                                                               Registered By                                 Description                                                      
======================= ================================================================================== ============================================= ================================================================ 
domainBoundaryIndicator integer_array                                                                                                                    (no description available)                                       
edgeList                geosx_InterObjectRelation< LvArray_ArrayOfSets< long, long, LvArray_ChaiBuffer > >                                               (no description available)                                       
elemList                LvArray_ArrayOfArrays< long, long, LvArray_ChaiBuffer >                                                                          (no description available)                                       
elemRegionList          LvArray_ArrayOfArrays< long, long, LvArray_ChaiBuffer >                                                                          (no description available)                                       
elemSubRegionList       LvArray_ArrayOfArrays< long, long, LvArray_ChaiBuffer >                                                                          (no description available)                                       
ghostRank               integer_array                                                                                                                    (no description available)                                       
globalToLocalMap        geosx_mapBase< long long, long, std_integral_constant< bool, false > >                                                           (no description available)                                       
isExternal              integer_array                                                                                                                    (no description available)                                       
localToGlobalMap        globalIndex_array                                                                                                                Array that contains a map from localIndex to globalIndex.        
parentEdgeGlobalIndex   globalIndex_array                                                                                                                (no description available)                                       
referencePosition       real64_array2d                                                                                                                   (no description available)                                       
parentEdgeIndex         localIndex_array                                                                   :ref:`DATASTRUCTURE_EmbeddedSurfaceGenerator` Index of parent edge within the mesh object it is registered on. 
neighborData            node                                                                                                                             :ref:`DATASTRUCTURE_neighborData`                                
sets                    node                                                                                                                             :ref:`DATASTRUCTURE_sets`                                        
======================= ================================================================================== ============================================= ================================================================ 


