

======================= =============================================================================== ============================================ ============================================================== 
Name                    Type                                                                            Registered By                                Description                                                    
======================= =============================================================================== ============================================ ============================================================== 
ReferencePosition       real64_array2d                                                                                                               (no description available)                                     
domainBoundaryIndicator integer_array                                                                                                                (no description available)                                     
edgeList                geos_InterObjectRelation< LvArray_ArrayOfSets< int, int, LvArray_ChaiBuffer > >                                              (no description available)                                     
elemList                LvArray_ArrayOfArrays< int, int, LvArray_ChaiBuffer >                                                                        (no description available)                                     
elemRegionList          LvArray_ArrayOfArrays< int, int, LvArray_ChaiBuffer >                                                                        (no description available)                                     
elemSubRegionList       LvArray_ArrayOfArrays< int, int, LvArray_ChaiBuffer >                                                                        (no description available)                                     
faceList                geos_InterObjectRelation< LvArray_ArrayOfSets< int, int, LvArray_ChaiBuffer > >                                              (no description available)                                     
ghostRank               integer_array                                                                                                                (no description available)                                     
globalToLocalMap        geos_mapBase< long long, int, std_integral_constant< bool, false > >                                                         (no description available)                                     
isExternal              integer_array                                                                                                                (no description available)                                     
localToGlobalMap        globalIndex_array                                                                                                            Array that contains a map from localIndex to globalIndex.      
primaryField            real64_array                                                                                                                 Primary field variable                                         
aSigma                  real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Product of direct effect parameter a and initial normal stress 
normalStress            real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Normal stress acting on the fault                              
normalStressRate        real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Normal stress rate acting on the fault                         
pressure                real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Pore pressure                                                  
pressureRate            real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Pore pressure rate                                             
shearStress             real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Shear stress acting on the fault                               
shearStressRate         real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Shear stress rate acting on the fault                          
t_a                     real64_array                                                                    :ref:`DATASTRUCTURE_DieterichSeismicityRate` Dieterich constitutive relaxation time of seismicity rate      
neighborData            node                                                                                                                         :ref:`DATASTRUCTURE_neighborData`                              
sets                    node                                                                                                                         :ref:`DATASTRUCTURE_sets`                                      
======================= =============================================================================== ============================================ ============================================================== 


