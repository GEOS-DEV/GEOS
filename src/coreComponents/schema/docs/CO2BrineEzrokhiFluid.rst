

==================== ============ ======== ============================================================================== 
Name                 Type         Default  Description                                                                    
==================== ============ ======== ============================================================================== 
componentMolarWeight real64_array {0}      Component molar weights                                                        
componentNames       string_array {}       List of component names                                                        
flashModelParaFile   path         required Name of the file defining the parameters of the flash model                    
name                 string       required A name is required for any non-unique nodes                                    
phaseNames           string_array {}       List of fluid phases                                                           
phasePVTParaFiles    path_array   required Names of the files defining the parameters of the viscosity and density models 
==================== ============ ======== ============================================================================== 


