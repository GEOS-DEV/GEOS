

=============== ======================================================== ======== =================================================================== 
Name            Type                                                     Default  Description                                                         
=============== ======================================================== ======== =================================================================== 
checkFrequency  integer                                                  10       MsRSB basis smoothing convergence check frequency                   
maxIter         integer                                                  100      Maximum number of MsRSB basis smoothing iterations                  
numLayers       integer                                                  3        Number of extra layers in support region (for supportType = layers) 
relaxation      real64                                                   0.666667 MsRSB basis smoothing iteration relaxation parameter                
supportType     geos_LinearSolverParameters_Multiscale_MsRSB_SupportType matching Type of support region construction algorithm                       
tolerance       real64                                                   0.001    MsRSB basis smoothing iteration tolerance                           
updateFrequency integer                                                  10       MsRSB basis smoothing coarse operator update frequency              
=============== ======================================================== ======== =================================================================== 


