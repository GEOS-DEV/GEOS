

======================================== ============ ======== ==================================================================================================== 
Name                                     Type         Default  Description                                                                                          
======================================== ============ ======== ==================================================================================================== 
bcApplicationTableName                   string                Name of table that specifies the on/off application of the bc.                                       
beginTime                                real64       -1e+99   time at which BC will start being applied.                                                           
componentFractionVsElevationTableNames   string_array {}       Names of the tables specifying the relationship (component fraction vs elevation) for each component 
componentNames                           string_array {}       Names of the fluid components                                                                        
datumElevation                           real64       required Datum elevation [m]                                                                                  
datumPressure                            real64       required Datum pressure [Pa]                                                                                  
direction                                R1Tensor     {0,0,0}  Direction to apply boundary condition to                                                             
endTime                                  real64       1e+99    time at which bc will stop being applied                                                             
equilibrationTolerance                   real64       0.001    Tolerance in the fixed-point iteration scheme used for hydrostatic initialization                    
functionName                             string                Name of function that specifies variation of the BC                                                  
initialPhaseName                         string                Name of the phase initially saturating the reservoir                                                 
maxNumberOfEquilibrationIterations       localIndex   5        Maximum number of equilibration iterations                                                           
name                                     string       required A name is required for any non-unique nodes                                                          
numberOfPointsInHydrostaticPressureTable localIndex   100      Number of points in the hydrostatic pressure table constructed internally                            
objectPath                               string                Path to the target field                                                                             
scale                                    real64       0        Scale factor for value of BC.                                                                        
setNames                                 string_array required Name of sets that boundary condition is applied to.                                                  
temperatureVsElevationTableName          string                Name of the table specifying the relationship (temperature [K] vs elevation)                         
======================================== ============ ======== ==================================================================================================== 


