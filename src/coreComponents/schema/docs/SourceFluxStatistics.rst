

============== ================== ======== =============================================================================================================================================== 
Name           Type               Default  Description                                                                                                                                     
============== ================== ======== =============================================================================================================================================== 
flowSolverName groupNameRef       required Name of the flow solver                                                                                                                         
fluxNames      groupNameRef_array required Name(s) array of the SourceFlux(s) for which we want the statistics. Use "all" to target all SourceFlux.                                        
logLevel       integer            0        | Log level                                                                                                                                       
                                           | - Log Level 1 outputs the sum of all SourceFlux(s) produced rate & mass,                                                                        
                                           | - Log Level 2 outputs detailed values for each SourceFlux.                                                                                      
name           groupName          required A name is required for any non-unique nodes                                                                                                     
============== ================== ======== =============================================================================================================================================== 


