

==================== ========= ============ ================================================================================================== 
Name                 Type      Default      Description                                                                                        
==================== ========= ============ ================================================================================================== 
dDrainedTEC_dT       real64    0            Derivative of the Thermal Expansion Coefficient of the Solid Rock Frame w.r.t. temperature [1/K^2] 
defaultBulkModulus   real64    -1           Default Bulk Modulus Parameter                                                                     
defaultDensity       real64    required     Default Material Density [Kg/cm^3]                                                                 
defaultDrainedTEC    real64    0            Default Linear Thermal Expansion Coefficient of the Solid Rock Frame [1/K]                         
defaultPoissonRatio  real64    -1           Default Poisson's Ratio                                                                            
defaultShearModulus  real64    -1           Default Shear Modulus Parameter                                                                    
defaultYieldStress   real64    1.79769e+308 Default yield stress                                                                               
defaultYoungModulus  real64    -1           Default Young's Modulus                                                                            
drainedTECTableName  string                 Name of the Thermal Expansion Coefficient table                                                    
name                 groupName required     A name is required for any non-unique nodes                                                        
referenceTemperature real64    0            Reference temperature at which the default Thermal Expansion Coefficient is defined [K]            
==================== ========= ============ ================================================================================================== 


