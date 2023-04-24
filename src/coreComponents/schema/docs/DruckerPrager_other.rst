

=========================== ===================================================================================== ========================================== 
Name                        Type                                                                                  Description                                
=========================== ===================================================================================== ========================================== 
bulkModulus                 real64_array                                                                          Elastic Bulk Modulus Field                 
cohesion                    real64_array2d                                                                        New cohesion state                         
density                     real64_array2d                                                                        Material Density                           
dilation                    real64_array                                                                          Plastic potential slope                    
friction                    real64_array                                                                          Yield surface slope                        
hardening                   real64_array                                                                          Hardening rate                             
oldCohesion                 real64_array2d                                                                        Old cohesion state                         
oldStress                   LvArray_Array< double, 3, camp_int_seq< long, 2l, 1l, 0l >, int, LvArray_ChaiBuffer > Previous Material Stress                   
shearModulus                real64_array                                                                          Elastic Shear Modulus Field                
stress                      LvArray_Array< double, 3, camp_int_seq< long, 2l, 1l, 0l >, int, LvArray_ChaiBuffer > Current Material Stress                    
thermalExpansionCoefficient real64_array                                                                          Linear Thermal Expansion Coefficient Field 
=========================== ===================================================================================== ========================================== 


