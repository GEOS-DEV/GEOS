

================ ======================== ======== ============================================================================= 
Name             Type                     Default  Description                                                                   
================ ======================== ======== ============================================================================= 
childDirectory   string                            Child directory path                                                          
format           geosx_vtk_VTKOutputMode  binary   Output data format.  Valid options: ``binary``, ``ascii``                     
name             string                   required A name is required for any non-unique nodes                                   
outputRegionType geosx_vtk_VTKRegionTypes all      Output region types.  Valid options: ``cell``, ``well``, ``surface``, ``all`` 
parallelThreads  integer                  1        Number of plot files.                                                         
plotFileRoot     string                   VTK      Name of the root file for this output.                                        
plotLevel        integer                  1        Level detail plot. Only fields with lower of equal plot level will be output. 
writeFEMFaces    integer                  0        (no description available)                                                    
================ ======================== ======== ============================================================================= 


