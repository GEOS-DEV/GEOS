

======================== ======= ====================================================================================================================== 
Name                     Type    Description                                                                                                            
======================== ======= ====================================================================================================================== 
beginFromRestart         integer Flag to indicate restart run.                                                                                          
inputFileName            string  Name of the input xml file.                                                                                            
outputDirectory          string  Directory in which to put the output files, if not specified defaults to the current directory.                        
overridePartitionNumbers integer Flag to indicate partition number override                                                                             
problemName              string  Used in writing the output files, if not specified defaults to the name of the input file.                             
restartFileName          string  Name of the restart file.                                                                                              
schemaFileName           string  Name of the output schema                                                                                              
suppressPinned           integer Whether to disallow using pinned memory allocations for MPI communication buffers.                                     
useNonblockingMPI        integer Whether to prefer using non-blocking MPI communication where implemented (results in non-deterministic DOF numbering). 
xPartitionsOverride      integer Number of partitions in the x-direction                                                                                
yPartitionsOverride      integer Number of partitions in the y-direction                                                                                
zPartitionsOverride      integer Number of partitions in the z-direction                                                                                
======================== ======= ====================================================================================================================== 


