

=========================== ===================================================================================== ========================================== 
Name                        Type                                                                                  Description                                
=========================== ===================================================================================== ========================================== 
bulkModulus                 real64_array                                                                          Elastic Bulk Modulus Field                 
density                     real64_array2d                                                                        Material Density                           
dilationRatio               real64_array                                                                          Plastic potential slope ratio              
hardening                   real64_array                                                                          Hardening parameter                        
initialFriction             real64_array                                                                          Initial yield surface slope                
oldStateVariable            real64_array2d                                                                        Old equivalent plastic shear strain        
oldStress                   LvArray_Array< double, 3, camp_int_seq< long, 2l, 1l, 0l >, int, LvArray_ChaiBuffer > Previous Material Stress                   
pressureIntercept           real64_array                                                                          Pressure point at cone vertex              
residualFriction            real64_array                                                                          Residual yield surface slope               
shearModulus                real64_array                                                                          Elastic Shear Modulus Field                
stateVariable               real64_array2d                                                                        New equivalent plastic shear strain        
stress                      LvArray_Array< double, 3, camp_int_seq< long, 2l, 1l, 0l >, int, LvArray_ChaiBuffer > Current Material Stress                    
thermalExpansionCoefficient real64_array                                                                          Linear Thermal Expansion Coefficient Field 
=========================== ===================================================================================== ========================================== 


