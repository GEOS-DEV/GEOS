

================= ====== ======== ================================================================================== 
Name              Type   Default  Description                                                                        
================= ====== ======== ================================================================================== 
fieldName         string required Name of primary solution field                                                     
boundaryFieldName string          Name of boundary (face) field                                                      
fractureRegions   string          Names of the fracture region that will have a fracture stencil generated for them. 
coefficientName   string required Name of coefficient field                                                          
areaRelTol        real64 1e-08    Relative tolerance for area calculations.                                          
name              string required A name is required for any non-unique nodes                                        
================= ====== ======== ================================================================================== 


