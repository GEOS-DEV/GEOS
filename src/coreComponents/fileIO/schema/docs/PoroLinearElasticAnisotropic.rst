

================= ============== ======== =============================================================== 
Name              Type           Default  Description                                                     
================= ============== ======== =============================================================== 
BiotCoefficient   real64         0        Biot's coefficient                                              
compressibility   real64         -1       Fluid Compressibilty                                            
defaultDensity    real64         required Default Material Density                                        
defaultStiffness  real64_array2d required Default Elastic Stiffness Tensor in Voigt notation (6x6 matrix) 
name              string         required A name is required for any non-unique nodes                     
referencePressure real64         0        ReferencePressure                                               
================= ============== ======== =============================================================== 


