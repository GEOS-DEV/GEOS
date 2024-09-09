Datastructure: ElementRegions
=============================

======================= ================================================================= ========================================================= 
Name                    Type                                                              Description                                               
======================= ================================================================= ========================================================= 
domainBoundaryIndicator integer_array                                                     (no description available)                                
ghostRank               integer_array                                                     (no description available)                                
globalToLocalMap        geos_mapBase<long long, int, std_integral_constant<bool, false> > (no description available)                                
isExternal              integer_array                                                     (no description available)                                
localMaxGlobalIndex     globalIndex                                                       (no description available)                                
localToGlobalMap        globalIndex_array                                                 Array that contains a map from localIndex to globalIndex. 
maxGlobalIndex          globalIndex                                                       (no description available)                                
CellElementRegion       node                                                              :ref:`DATASTRUCTURE_CellElementRegion`                    
SurfaceElementRegion    node                                                              :ref:`DATASTRUCTURE_SurfaceElementRegion`                 
WellElementRegion       node                                                              :ref:`DATASTRUCTURE_WellElementRegion`                    
elementRegionsGroup     node                                                              :ref:`DATASTRUCTURE_elementRegionsGroup`                  
neighborData            node                                                              :ref:`DATASTRUCTURE_neighborData`                         
sets                    node                                                              :ref:`DATASTRUCTURE_sets`                                 
======================= ================================================================= ========================================================= 


