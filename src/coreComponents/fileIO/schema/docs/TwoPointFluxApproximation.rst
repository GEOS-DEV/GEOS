

================= ============ ======== =========================================== 
Name              Type         Default  Description                                 
================= ============ ======== =========================================== 
fieldName         string       required Name of primary solution field              
boundaryFieldName string                Name of boundary (face) field               
coefficientName   string       required Name of coefficient field                   
targetRegions     string_array          List of regions to build the stencil for    
areaRelTol        real64       1e-08    Relative tolerance for area calculations.   
name              string       required A name is required for any non-unique nodes 
================= ============ ======== =========================================== 


