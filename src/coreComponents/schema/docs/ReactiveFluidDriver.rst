

================== ============ ======== ============================================================================================================================================================= 
Name               Type         Default  Description                                                                                                                                                   
================== ============ ======== ============================================================================================================================================================= 
baseline           path         none     Baseline file                                                                                                                                                 
feedComposition    real64_array required Feed composition array: total concentration of the primary species                                                                                            
fluid              groupNameRef required Fluid to test                                                                                                                                                 
logLevel           integer      0        | Sets the level of information to write in the standard output (the console typically).                                                                        
                                         | A level of 0 outputs minimal information, higher levels require more.                                                                                         
name               groupName    required A name is required for any non-unique nodes                                                                                                                   
output             string       none     Output file                                                                                                                                                   
pressureControl    groupNameRef required Function controlling pressure time history                                                                                                                    
steps              integer      required Number of load steps to take                                                                                                                                  
temperatureControl groupNameRef required Function controlling temperature time history                                                                                                                 
================== ============ ======== ============================================================================================================================================================= 


