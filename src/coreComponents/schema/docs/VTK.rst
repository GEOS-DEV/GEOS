

=========================== ======================= ======== ======================================================================================================================================================================================== 
Name                        Type                    Default  Description                                                                                                                                                                              
=========================== ======================= ======== ======================================================================================================================================================================================== 
childDirectory              string                           Child directory path                                                                                                                                                                     
fieldNames                  string_array            {}       Names of the fields to output. If this attribute is specified, GEOSX outputs all the fields specified by the user, regardless of their `plotLevel`                                       
format                      geos_vtk_VTKOutputMode  binary   Output data format.  Valid options: ``binary``, ``ascii``                                                                                                                                
name                        string                  required A name is required for any non-unique nodes                                                                                                                                              
onlyPlotSpecifiedFieldNames integer                 0        If this flag is equal to 1, then we only plot the fields listed in `fieldNames`. Otherwise, we plot all the fields with the required `plotLevel`, plus the fields listed in `fieldNames` 
outputRegionType            geos_vtk_VTKRegionTypes all      Output region types.  Valid options: ``cell``, ``well``, ``surface``, ``all``                                                                                                            
parallelThreads             integer                 1        Number of plot files.                                                                                                                                                                    
plotFileRoot                string                  VTK      Name of the root file for this output.                                                                                                                                                   
plotLevel                   integer                 1        Level detail plot. Only fields with lower of equal plot level will be output.                                                                                                            
writeFEMFaces               integer                 0        (no description available)                                                                                                                                                               
writeGhostCells             integer                 0        Should the vtk files contain the ghost cells or not.                                                                                                                                     
=========================== ======================= ======== ======================================================================================================================================================================================== 


