

===================== ==================================================================================================================================================== ======== ================================================================================== 
Name                  Type                                                                                                                                                 Default  Description                                                                        
===================== ==================================================================================================================================================== ======== ================================================================================== 
areaRelTol            real64                                                                                                                                               1e-08    Relative tolerance for area calculations.                                          
coefficientModelNames geosx_mapBase< std_string, LvArray_Array< std_string, 1, camp_int_seq< long, 0l >, long, LvArray_ChaiBuffer >, std_integral_constant< bool, true > >          List of constitutive models that contain the coefficient used to build the stencil 
coefficientName       string                                                                                                                                               required Name of coefficient field                                                          
fieldName             string                                                                                                                                               required Name of primary solution field                                                     
meanPermCoefficient   real64                                                                                                                                               1        (no description available)                                                         
name                  string                                                                                                                                               required A name is required for any non-unique nodes                                        
targetRegions         geosx_mapBase< std_string, LvArray_Array< std_string, 1, camp_int_seq< long, 0l >, long, LvArray_ChaiBuffer >, std_integral_constant< bool, true > >          List of regions to build the stencil for                                           
usePEDFM              integer                                                                                                                                              0        (no description available)                                                         
===================== ==================================================================================================================================================== ======== ================================================================================== 


