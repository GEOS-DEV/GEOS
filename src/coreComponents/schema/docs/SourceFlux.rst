

====================== ================== ======== ============================================================================================================================================================= 
Name                   Type               Default  Description                                                                                                                                                   
====================== ================== ======== ============================================================================================================================================================= 
bcApplicationTableName groupNameRef                Name of table that specifies the on/off application of the boundary condition.                                                                                
beginTime              real64             -1e+99   Time at which the boundary condition will start being applied.                                                                                                
component              integer            -1       Component of field (if tensor) to apply boundary condition to.                                                                                                
direction              R1Tensor           {0,0,0}  Direction to apply boundary condition to.                                                                                                                     
endTime                real64             1e+99    Time at which the boundary condition will stop being applied.                                                                                                 
functionName           groupNameRef                Name of function that specifies variation of the boundary condition.                                                                                          
initialCondition       integer            0        Boundary condition is applied as an initial condition.                                                                                                        
logLevel               integer            0        | Sets the level of information to write in the standard output (the console typically).                                                                        
                                                   | A level of 0 outputs minimal information, higher levels require more.                                                                                         
name                   groupName          required A name is required for any non-unique nodes                                                                                                                   
objectPath             groupNameRef                Path to the target field                                                                                                                                      
scale                  real64             0        Scale factor for value of the boundary condition.                                                                                                             
setNames               groupNameRef_array required Name of sets that boundary condition is applied to.                                                                                                           
====================== ================== ======== ============================================================================================================================================================= 


