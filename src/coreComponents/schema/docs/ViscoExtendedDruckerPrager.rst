

============================ ========= ======== ================================================================================================== 
Name                         Type      Default  Description                                                                                        
============================ ========= ======== ================================================================================================== 
dDrainedTEC_dT               real64    0        Derivative of the Thermal Expansion Coefficient of the Solid Rock Frame w.r.t. temperature [1/K^2] 
defaultBulkModulus           real64    -1       Default Bulk Modulus Parameter                                                                     
defaultCohesion              real64    0        Initial cohesion                                                                                   
defaultDensity               real64    required Default Material Density [Kg/cm^3]                                                                 
defaultDilationRatio         real64    1        Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle)                               
defaultDrainedTEC            real64    0        Default Linear Thermal Expansion Coefficient of the Solid Rock Frame [1/K]                         
defaultHardening             real64    0        Hardening parameter (hardening rate is faster for smaller values)                                  
defaultInitialFrictionAngle  real64    30       Initial friction angle (degrees)                                                                   
defaultPoissonRatio          real64    -1       Default Poisson's Ratio                                                                            
defaultResidualFrictionAngle real64    30       Residual friction angle (degrees)                                                                  
defaultShearModulus          real64    -1       Default Shear Modulus Parameter                                                                    
defaultYoungModulus          real64    -1       Default Young's Modulus                                                                            
drainedTECTableName          string             Name of the Thermal Expansion Coefficient table                                                    
name                         groupName required A name is required for any non-unique nodes                                                        
referenceTemperature         real64    0        Reference temperature at which the default Thermal Expansion Coefficient is defined [K]            
relaxationTime               real64    required Relaxation time                                                                                    
============================ ========= ======== ================================================================================================== 


