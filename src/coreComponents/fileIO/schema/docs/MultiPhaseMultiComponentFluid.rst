

==================== ============ ======== ================================================================= 
Name                 Type         Default  Description                                                       
==================== ============ ======== ================================================================= 
componentMolarWeight real64_array {0}      Component molar weights                                           
componentNames       string_array {}       List of component names                                           
flashModelParaFile   path         required name of the filen including flash calculation function parameters 
name                 string       required A name is required for any non-unique nodes                       
phaseNames           string_array {}       List of fluid phases                                              
phasePVTParaFiles    path_array   required List of the names of the files including PVT function parameters  
==================== ============ ======== ================================================================= 


