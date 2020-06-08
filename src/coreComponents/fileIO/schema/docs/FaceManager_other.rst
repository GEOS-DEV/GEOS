

======================= ====================================================================== ========================================= ================================================================================== 
Name                    Type                                                                   Registered By                             Description                                                                        
======================= ====================================================================== ========================================= ================================================================================== 
K_IC                    r1_array                                                                                                         K_IC on the face                                                                   
SIFonFace               real64_array                                                                                                     SIF on the face                                                                    
boundaryFaceDensity     real64_array2d                                                                                                   (no description available)                                                         
boundaryFaceMobility    real64_array                                                                                                     (no description available)                                                         
boundaryFacePressure    real64_array                                                                                                     (no description available)                                                         
boundaryFaceViscosity   real64_array2d                                                                                                   (no description available)                                                         
childIndex              localIndex_array                                                                                                 child index of the face.                                                           
degreeFromCrackTip      integer_array                                                                                                    degree of connectivity separation from crack tip.                                  
domainBoundaryIndicator integer_array                                                                                                    (no description available)                                                         
edgeList                geosx_InterObjectRelation< LvArray_ArrayOfArrays< long, long > >                                                 (no description available)                                                         
elemList                localIndex_array2d                                                                                               (no description available)                                                         
elemRegionList          localIndex_array2d                                                                                               (no description available)                                                         
elemSubRegionList       localIndex_array2d                                                                                               (no description available)                                                         
faceArea                real64_array                                                                                                     (no description available)                                                         
faceCenter              real64_array2d                                                                                                   (no description available)                                                         
faceNormal              real64_array2d                                                                                                   (no description available)                                                         
faceRotationMatrix      real64_array3d                                                                                                   (no description available)                                                         
ghostRank               integer_array                                                                                                    (no description available)                                                         
globalToLocalMap        geosx_mapBase< long long, long, std_integral_constant< bool, false > >                                           (no description available)                                                         
gravityCoefficient      real64_array                                                                                                     (no description available)                                                         
isExternal              integer_array                                                                                                    (no description available)                                                         
isFaceSeparable         integer_array                                                                                                    A flag to mark if the face is separable                                            
localToGlobalMap        globalIndex_array                                                                                                Array that contains a map from localIndex to globalIndex.                          
nodeList                geosx_InterObjectRelation< LvArray_ArrayOfArrays< long, long > >                                                 (no description available)                                                         
parentIndex             localIndex_array                                                                                                 Parent index of the face.                                                          
primaryCandidateFace    localIndex_array                                                                                                 The face that has the highest score for splitability                               
ruptureState            integer_array                                                                                                    Rupture state of the face.0=not ready for rupture. 1=ready for rupture. 2=ruptured 
ruptureTime             real64_array                                                                                                     Time that the face was ruptured.                                                   
deltaFacePressure       real64_array                                                           :ref:`DATASTRUCTURE_SinglePhaseHybridFVM` An array that holds the accumulated pressure updates at the faces.                 
facePressure            real64_array                                                           :ref:`DATASTRUCTURE_SinglePhaseHybridFVM` An array that holds the pressures at the faces.                                    
neighborData            node                                                                                                             :ref:`DATASTRUCTURE_neighborData`                                                  
sets                    node                                                                                                             :ref:`DATASTRUCTURE_sets`                                                          
======================= ====================================================================== ========================================= ================================================================================== 


