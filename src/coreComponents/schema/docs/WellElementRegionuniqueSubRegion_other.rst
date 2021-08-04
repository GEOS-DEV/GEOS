

========================== ============================================================================================================= ========================================================= 
Name                       Type                                                                                                          Description                                               
========================== ============================================================================================================= ========================================================= 
domainBoundaryIndicator    integer_array                                                                                                 (no description available)                                
elementCenter              real64_array2d                                                                                                (no description available)                                
elementVolume              real64_array                                                                                                  (no description available)                                
ghostRank                  integer_array                                                                                                 (no description available)                                
globalToLocalMap           geosx_mapBase< long long, long, std_integral_constant< bool, false > >                                        (no description available)                                
isExternal                 integer_array                                                                                                 (no description available)                                
localToGlobalMap           globalIndex_array                                                                                             Array that contains a map from localIndex to globalIndex. 
nextWellElementIndex       localIndex_array                                                                                              (no description available)                                
nextWellElementIndexGlobal localIndex_array                                                                                              (no description available)                                
nodeList                   geosx_InterObjectRelation< LvArray_Array< long, 2, camp_int_seq< long, 0l, 1l >, long, LvArray_ChaiBuffer > > (no description available)                                
numEdgesPerElement         localIndex                                                                                                    (no description available)                                
numFacesPerElement         localIndex                                                                                                    (no description available)                                
numNodesPerElement         localIndex                                                                                                    (no description available)                                
radius                     real64_array                                                                                                  (no description available)                                
topRank                    integer                                                                                                       (no description available)                                
topWellElementIndex        localIndex                                                                                                    (no description available)                                
wellControlsName           string                                                                                                        (no description available)                                
ConstitutiveModels         node                                                                                                          :ref:`DATASTRUCTURE_ConstitutiveModels`                   
neighborData               node                                                                                                          :ref:`DATASTRUCTURE_neighborData`                         
sets                       node                                                                                                          :ref:`DATASTRUCTURE_sets`                                 
wellElementSubRegion       node                                                                                                          :ref:`DATASTRUCTURE_wellElementSubRegion`                 
========================== ============================================================================================================= ========================================================= 


