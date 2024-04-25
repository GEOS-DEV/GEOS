

=========================== ============ ======== ============================================================================================================================================================= 
Name                        Type         Default  Description                                                                                                                                                   
=========================== ============ ======== ============================================================================================================================================================= 
logLevel                    integer      1        | Sets the level of information to write in the standard output (the console typically).                                                                        
                                                  | A level of 0 outputs minimal information, higher levels require more.                                                                                         
name                        groupName    required A name is required for any non-unique nodes                                                                                                                   
performStressInitialization integer      required Flag to indicate that the solver is going to perform stress initialization                                                                                    
poromechanicsSolverName     groupNameRef required Name of the poromechanics solver                                                                                                                              
=========================== ============ ======== ============================================================================================================================================================= 


