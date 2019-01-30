

====================== ============ ======= ============================================================== 
Name                   Type         Default Description                                                    
====================== ============ ======= ============================================================== 
setNames               string_array         Name of sets that boundary condition is applied to.            
objectPath             string               Path to the target field                                       
fieldName              string               Name of field that boundary condition is applied to.           
component              integer      0       Component of field (if tensor) to apply boundary condition to  
direction              R1Tensor     0 0 0   Direction to apply boundary condition to                       
functionName           string               Name of function that specifies variation of the BC            
bcApplicationTableName string               Name of table that specifies the on/off application of the bc. 
scale                  real64       0       Scale factor for value of BC.                                  
initialCondition       integer      0       BC is applied as an initial condition.                         
beginTime              real64       -1e+99  time at which BC will start being applied.                     
endTime                real64       1e+99   time at which bc will stop being applied                       
name                   string               A name is required for any non-unique nodes                    
====================== ============ ======= ============================================================== 


