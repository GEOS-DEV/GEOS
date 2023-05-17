

============== ============= ======== ======================================================================================================================================================= 
Name           Type          Default  Description                                                                                                                                             
============== ============= ======== ======================================================================================================================================================= 
args           string                 Any extra command line arguments to pass to GEOSX.                                                                                                      
autoPartition  string                 May be 'Off' or 'On', if 'On' partitioning arguments are created automatically. Default is Off.                                                         
meshSizes      integer_array {0}      The target number of elements in the internal mesh (per-process for weak scaling, globally for strong scaling) default doesn't modify the internalMesh. 
name           string        required The name of this benchmark.                                                                                                                             
nodes          integer       0        The number of nodes needed to run the base benchmark, default is 1.                                                                                     
scaleList      integer_array {0}      The scales at which to run the problem ( scale * nodes * tasksPerNode ).                                                                                
scaling        string                 Whether to run a scaling, and which type of scaling to run.                                                                                             
tasksPerNode   integer       required The number of tasks per node to run the benchmark with.                                                                                                 
threadsPerTask integer       0        The number of threads per task to run the benchmark with.                                                                                               
timeLimit      integer       0        The time limit of the benchmark.                                                                                                                        
============== ============= ======== ======================================================================================================================================================= 


