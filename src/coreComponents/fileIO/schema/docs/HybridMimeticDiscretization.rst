

================ ============ ======== =================================================== 
Name             Type         Default  Description                                         
================ ============ ======== =================================================== 
areaRelTol       real64       1e-08    Relative tolerance for area calculations.           
coefficientName  string       required Name of coefficient field                           
fieldName        string       required Name of primary solution field                      
innerProductType string       required Type of inner product used in the hybrid FVM solver 
name             string       required A name is required for any non-unique nodes         
targetRegions    string_array {}       List of regions to build the stencil for            
================ ============ ======== =================================================== 


