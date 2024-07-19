

====================== ========= ======== ================================================================================================== 
Name                   Type      Default  Description                                                                                        
====================== ========= ======== ================================================================================================== 
compressiveStrength    real64    0        Compressive strength from the uniaxial compression test                                            
criticalFractureEnergy real64    required Critical fracture energy                                                                           
criticalStrainEnergy   real64    required Critical stress in a 1d tension test                                                               
dDrainedTEC_dT         real64    0        Derivative of the Thermal Expansion Coefficient of the Solid Rock Frame w.r.t. temperature [1/K^2] 
defaultBulkModulus     real64    -1       Default Bulk Modulus Parameter                                                                     
defaultDensity         real64    required Default Material Density [Kg/cm^3]                                                                 
defaultDrainedTEC      real64    0        Default Linear Thermal Expansion Coefficient of the Solid Rock Frame [1/K]                         
defaultPoissonRatio    real64    -1       Default Poisson's Ratio                                                                            
defaultShearModulus    real64    -1       Default Shear Modulus Parameter                                                                    
defaultYoungModulus    real64    -1       Default Young's Modulus                                                                            
degradationLowerLimit  real64    0        The lower limit of the degradation function                                                        
deltaCoefficient       real64    -1       Coefficient in the calculation of the external driving force                                       
drainedTECTableName    string             Name of the Thermal Expansion Coefficient table                                                    
extDrivingForceFlag    integer   0        Whether to have external driving force. Can be 0 or 1                                              
lengthScale            real64    required Length scale l in the phase-field equation                                                         
name                   groupName required A name is required for any non-unique nodes                                                        
referenceTemperature   real64    0        Reference temperature at which the default Thermal Expansion Coefficient is defined [K]            
tensileStrength        real64    0        Tensile strength from the uniaxial tension test                                                    
====================== ========= ======== ================================================================================================== 


