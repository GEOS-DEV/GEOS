

============================ ================== =============== ========================================================================================================================= 
Name                         Type               Default         Description                                                                                                               
============================ ================== =============== ========================================================================================================================= 
checkPVTTablesRanges         integer            1               Enable (1) or disable (0) an error when the input pressure or temperature of the PVT tables is out of range.              
componentAcentricFactor      real64_array       required        Component acentric factors                                                                                                
componentBinaryCoeff         real64_array2d     {{0}}           Table of binary interaction coefficients                                                                                  
componentCriticalPressure    real64_array       required        Component critical pressures                                                                                              
componentCriticalTemperature real64_array       required        Component critical temperatures                                                                                           
componentCriticalVolume      real64_array       {0}             Component critical volumes                                                                                                
componentMolarWeight         real64_array       required        Component molar weights                                                                                                   
componentNames               string_array       required        List of component names                                                                                                   
componentVolumeShift         real64_array       {0}             Component volume shifts                                                                                                   
equationsOfState             string_array       required        | List of equation of state types for each phase. Valid options:                                                            
                                                                | * pr                                                                                                                      
                                                                | * srk                                                                                                                     
name                         groupName          required        A name is required for any non-unique nodes                                                                               
phaseNames                   groupNameRef_array required        List of fluid phases                                                                                                      
viscosityMixingRule          string             HerningZipperer | Viscosity mixing rule to be used for Lohrenz-Bray-Clark computation. Valid options:                                       
                                                                | * HerningZipperer                                                                                                         
                                                                | * Wilke                                                                                                                   
                                                                | * Brokaw                                                                                                                  
============================ ================== =============== ========================================================================================================================= 


