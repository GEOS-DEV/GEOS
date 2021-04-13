

============== ============= ======== =============================================================================================== 
Name           Type          Default  Description                                                                                     
============== ============= ======== =============================================================================================== 
args           string                 Any extra command line arguments to pass to GEOSX.                                              
autoPartition  string                 May be 'Off' or 'On', if 'On' partitioning arguments are created automatically. Default is Off. 
name           string        required The name of this benchmark.                                                                     
nodes          integer       required The number of nodes needed to run the benchmark.                                                
strongScaling  integer_array {0}      Repeat the benchmark N times, scaling the number of nodes in the benchmark by these values.     
tasksPerNode   integer       required The number of tasks per node to run the benchmark with.                                         
threadsPerTask integer       0        The number of threads per task to run the benchmark with.                                       
timeLimit      integer       0        The time limit of the benchmark.                                                                
============== ============= ======== =============================================================================================== 


