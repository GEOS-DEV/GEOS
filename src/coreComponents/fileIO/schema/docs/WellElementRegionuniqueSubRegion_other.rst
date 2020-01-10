

========================== ===================================================================================== ========================================================= 
Name                       Type                                                                                  Description                                               
========================== ===================================================================================== ========================================================= 
domainBoundaryIndicator    integer_array                                                                         (no description available)                                
elementCenter              r1_array                                                                              (no description available)                                
elementVolume              real64_array                                                                          (no description available)                                
ghostRank                  integer_array                                                                         (no description available)                                
globalToLocalMap           mapBase< long long, long, integral_constant< bool, false > >                          (no description available)                                
isExternal                 integer_array                                                                         (no description available)                                
localToGlobalMap           globalIndex_array                                                                     Array that contains a map from localIndex to globalIndex. 
nextWellElementIndex       localIndex_array                                                                      (no description available)                                
nextWellElementIndexGlobal localIndex_array                                                                      (no description available)                                
nodeList                   InterObjectRelation< Array< long, 2, int_seq< long, 0l, 1l >, long, NewChaiBuffer > > (no description available)                                
topRank                    integer                                                                               (no description available)                                
topWellElementIndex        localIndex                                                                            (no description available)                                
wellControlsName           string                                                                                (no description available)                                
ConstitutiveModels         node                                                                                  :ref:`DATASTRUCTURE_ConstitutiveModels`                   
neighborData               node                                                                                  :ref:`DATASTRUCTURE_neighborData`                         
sets                       node                                                                                  :ref:`DATASTRUCTURE_sets`                                 
wellElementSubRegion       node                                                                                  :ref:`DATASTRUCTURE_wellElementSubRegion`                 
========================== ===================================================================================== ========================================================= 


