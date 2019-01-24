
Element: CompositionalMultiphaseFluid
=====================================

============================ ============== ======= ======== ============================================== 
Name                         Type           Default Use      Description                                    
============================ ============== ======= ======== ============================================== 
componentNames               string_array           required List of component names                        
componentMolarWeight         real64_array           required Component molar weights                        
phaseNames                   string_array           required List of fluid phases                           
equationsOfState             string_array           required List of equation of state types for each phase 
componentCriticalPressure    real64_array           required Component critical pressures                   
componentCriticalTemperature real64_array           required Component critical temperatures                
componentAcentricFactor      real64_array           required Component acentric factors                     
componentVolumeShift         real64_array   0                Component volume shifts                        
componentBinaryCoeff         real64_array2d 0                Table of binary interaction coefficients       
name                         string                 required                                                
============================ ============== ======= ======== ============================================== 


