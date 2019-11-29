

======================= =============================================================== ================================================================================== 
Name                    Type                                                            Description                                                                        
======================= =============================================================== ================================================================================== 
K_IC                    r1_array                                                        K_IC on the face                                                                   
SIFonFace               real64_array                                                    SIF on the face                                                                    
childIndex              localIndex_array                                                child index of the face.                                                           
domainBoundaryIndicator integer_array                                                   (no description available)                                                         
edgeList                InterObjectRelation< ArrayOfArrays< long, long > >              (no description available)                                                         
elemList                localIndex_array2d                                              (no description available)                                                         
elemRegionList          localIndex_array2d                                              (no description available)                                                         
elemSubRegionList       localIndex_array2d                                              (no description available)                                                         
faceArea                real64_array                                                    (no description available)                                                         
faceCenter              r1_array                                                        (no description available)                                                         
faceDensity             real64_array2d                                                  (no description available)                                                         
faceMobility            real64_array                                                    (no description available)                                                         
faceNormal              r1_array                                                        (no description available)                                                         
facePressure            real64_array                                                    (no description available)                                                         
faceViscosity           real64_array2d                                                  (no description available)                                                         
ghostRank               integer_array                                                   (no description available)                                                         
globalToLocalMap        mapBase< long long, long, __1integral_constant< bool, false > > (no description available)                                                         
gravityDepth            real64_array                                                    (no description available)                                                         
isExternal              integer_array                                                   (no description available)                                                         
isFaceSeparable         integer_array                                                   A flag to mark if the face is separable                                            
localToGlobalMap        globalIndex_array                                               Array that contains a map from localIndex to globalIndex.                          
nodeList                InterObjectRelation< ArrayOfArrays< long, long > >              (no description available)                                                         
parentIndex             localIndex_array                                                Parent index of the face.                                                          
primaryCandidateFace    localIndex_array                                                The face that has the highest score for splitability                               
ruptureState            integer_array                                                   Rupture state of the face.0=not ready for rupture. 1=ready for rupture. 2=ruptured 
neighborData            node                                                            :ref:`DATASTRUCTURE_neighborData`                                                  
sets                    node                                                            :ref:`DATASTRUCTURE_sets`                                                          
======================= =============================================================== ================================================================================== 


