

==================== ================== ======== ============================================================================================================ 
Name                 Type               Default  Description                                                                                                  
==================== ================== ======== ============================================================================================================ 
checkPVTTablesRanges integer            1        Enable (1) or disable (0) an error when the input pressure or temperature of the PVT tables is out of range. 
componentMolarWeight real64_array       {0}      Component molar weights                                                                                      
componentNames       string_array       {}       List of component names                                                                                      
logLevel             integer            0        Log level                                                                                                    
name                 groupName          required A name is required for any non-unique nodes                                                                  
phaseNames           groupNameRef_array {}       List of fluid phases                                                                                         
==================== ================== ======== ============================================================================================================ 


