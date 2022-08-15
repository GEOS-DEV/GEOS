

====================== ============ ======== ============================================================ 
Name                   Type         Default  Description                                                  
====================== ============ ======== ============================================================ 
biotCoefficient        real64_array {0}      Biot coefficient                                             
compressiveStrength    real64       0        Compressive strength from the uniaxial compression test      
criticalFractureEnergy real64       required Critical fracture energy                                     
criticalStrainEnergy   real64       required Critical stress in a 1d tension test                         
defaultBulkModulus     real64       -1       Default Bulk Modulus Parameter                               
defaultDensity         real64       required Default Material Density                                     
defaultPoissonRatio    real64       -1       Default Poisson's Ratio                                      
defaultShearModulus    real64       -1       Default Shear Modulus Parameter                              
defaultYoungModulus    real64       -1       Default Young's Modulus                                      
deltaCoefficient       real64       0        Coefficient in the calculation of the external driving force 
extDrivingForceSwitch  string       required Whether to have external driving force. Can be True or False 
lengthScale            real64       required Length scale l in the phase-field equation                   
name                   string       required A name is required for any non-unique nodes                  
tensileStrength        real64       0        Tensile strength from the uniaxial tension test              
====================== ============ ======== ============================================================ 


