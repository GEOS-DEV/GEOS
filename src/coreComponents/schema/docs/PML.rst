

====================== ============ =========================================== ================================================================================================== 
Name                   Type         Default                                     Description                                                                                        
====================== ============ =========================================== ================================================================================================== 
bcApplicationTableName string                                                   Name of table that specifies the on/off application of the boundary condition.                     
beginTime              real64       -1e+99                                      Time at which the boundary condition will start being applied.                                     
component              integer      -1                                          Component of field (if tensor) to apply boundary condition to.                                     
direction              R1Tensor     {0,0,0}                                     Direction to apply boundary condition to.                                                          
endTime                real64       1e+99                                       Time at which the boundary condition will stop being applied.                                      
functionName           string                                                   Name of function that specifies variation of the boundary condition.                               
logLevel               integer      0                                           Log level                                                                                          
name                   string       required                                    A name is required for any non-unique nodes                                                        
objectPath             string                                                   Path to the target field                                                                           
reflectivity           real64       0.001                                       Desired reflectivity of the PML region, used to compute the damping profile                        
scale                  real64       0                                           Scale factor for value of the boundary condition.                                                  
setNames               string_array required                                    Name of sets that boundary condition is applied to.                                                
thicknessMaxXYZ        R1Tensor32b  {-1,-1,-1}                                  Thickness of the PML region, at right, back, and bottom sides, used to compute the damping profile 
thicknessMinXYZ        R1Tensor32b  {-1,-1,-1}                                  Thickness of the PML region, at left, front, and top sides, used to compute the damping profile    
waveSpeedMaxXYZ        R1Tensor32b  {-1,-1,-1}                                  Wave speed in the PML, at right, back, and bottom sides, used to compute the damping profile       
waveSpeedMinXYZ        R1Tensor32b  {-1,-1,-1}                                  Wave speed in the PML, at left, front, and top sides, used to compute the damping profile          
xMax                   R1Tensor32b  {1.79769e+308,1.79769e+308,1.79769e+308}    Maximum (x,y,z) coordinates of the inner PML boundaries                                            
xMin                   R1Tensor32b  {-1.79769e+308,-1.79769e+308,-1.79769e+308} Minimum (x,y,z) coordinates of the inner PML boundaries                                            
====================== ============ =========================================== ================================================================================================== 


