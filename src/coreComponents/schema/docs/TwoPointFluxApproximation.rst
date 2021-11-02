

===================== ============ ======== ================================================================================== 
Name                  Type         Default  Description                                                                        
===================== ============ ======== ================================================================================== 
areaRelTol            real64       1e-08    Relative tolerance for area calculations.                                          
coefficientModelNames string_array {}       List of constitutive models that contain the coefficient used to build the stencil 
coefficientName       string       required Name of coefficient field                                                          
fieldName             string       required Name of primary solution field                                                     
meanPermCoefficient   real64       1        (no description available)                                                         
name                  string       required A name is required for any non-unique nodes                                        
targetRegions         string_array {}       List of regions to build the stencil for                                           
usePEDFM              integer      0        (no description available)                                                         
===================== ============ ======== ================================================================================== 


