

=========================== ============ ======== ======================================================================================================================================================================================== 
Name                        Type         Default  Description                                                                                                                                                                              
=========================== ============ ======== ======================================================================================================================================================================================== 
childDirectory              string                Child directory path                                                                                                                                                                     
fieldNames                  string_array {}       Names of the fields to output. If this attribute is specified, GEOSX outputs all the fields specified by the user, regardless of their `plotLevel`                                       
name                        string       required A name is required for any non-unique nodes                                                                                                                                              
objectName                  string                The name of the object from which to retrieve field values.                                                                                                                              
onlyPlotSpecifiedFieldNames integer      0        If this flag is equal to 1, then we only plot the fields listed in `fieldNames`. Otherwise, we plot all the fields with the required `plotLevel`, plus the fields listed in `fieldNames` 
parallelThreads             integer      1        Number of plot files.                                                                                                                                                                    
parentMeshName              string                The name of the grid.                                                                                                                                                                    
parentMeshUUID              string                The UUID of the grid.                                                                                                                                                                    
plotFileName                string       RESQML   Name of the file for this output.                                                                                                                                                        
plotLevel                   integer      1        Level detail plot. Only fields with lower of equal plot level will be output.                                                                                                            
=========================== ============ ======== ======================================================================================================================================================================================== 


