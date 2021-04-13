

========================= ============ ======== ============================================================================================================================================================ 
Name                      Type         Default  Description                                                                                                                                                  
========================= ============ ======== ============================================================================================================================================================ 
gasOilRelPermTableNames   string_array {}       | List of relative permeability tables for the pair (gas phase, oil phase)                                                                                     
                                                | The expected format is "{ gasPermTableName, oilPermTableName }", in that order                                                                               
name                      string       required A name is required for any non-unique nodes                                                                                                                  
phaseNames                string_array required List of fluid phases                                                                                                                                         
waterOilRelPermTableNames string_array {}       | List of relative permeability tables for the pair (water phase, oil phase)                                                                                   
                                                | The expected format is "{ waterPermTableName, oilPermTableName }", in that order                                                                             
========================= ============ ======== ============================================================================================================================================================ 


