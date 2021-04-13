

============================ ============== ======== ============================================== 
Name                         Type           Default  Description                                    
============================ ============== ======== ============================================== 
componentAcentricFactor      real64_array   required Component acentric factors                     
componentBinaryCoeff         real64_array2d {{0}}    Table of binary interaction coefficients       
componentCriticalPressure    real64_array   required Component critical pressures                   
componentCriticalTemperature real64_array   required Component critical temperatures                
componentMolarWeight         real64_array   required Component molar weights                        
componentNames               string_array   required List of component names                        
componentVolumeShift         real64_array   {0}      Component volume shifts                        
equationsOfState             string_array   required List of equation of state types for each phase 
name                         string         required A name is required for any non-unique nodes    
phaseNames                   string_array   required List of fluid phases                           
============================ ============== ======== ============================================== 


