

=============================== ========================================================================================= ========================================================================================================================== 
Name                            Type                                                                                      Description                                                                                                                
=============================== ========================================================================================= ========================================================================================================================== 
dPhaseRelPerm_dPhaseVolFraction LvArray_Array< double, 4, camp_int_seq< long, 1l, 2l, 3l, 0l >, int, LvArray_ChaiBuffer > Derivative of phase relative permeability with respect to phase volume fraction                                            
phaseOrder                      integer_array                                                                             (no description available)                                                                                                 
phaseRelPerm                    LvArray_Array< double, 3, camp_int_seq< long, 1l, 2l, 0l >, int, LvArray_ChaiBuffer >     Phase relative permeability                                                                                                
phaseRelPerm_n                  LvArray_Array< double, 3, camp_int_seq< long, 1l, 2l, 0l >, int, LvArray_ChaiBuffer >     Phase relative permeability at previous time                                                                               
phaseTrappedVolFraction         LvArray_Array< double, 3, camp_int_seq< long, 1l, 2l, 0l >, int, LvArray_ChaiBuffer >     Phase trapped volume fraction                                                                                              
phaseTypes                      integer_array                                                                             (no description available)                                                                                                 
volFracScale                    real64                                                                                    Factor used to scale the phase relative permeability, defined as: one minus the sum of the phase minimum volume fractions. 
=============================== ========================================================================================= ========================================================================================================================== 


