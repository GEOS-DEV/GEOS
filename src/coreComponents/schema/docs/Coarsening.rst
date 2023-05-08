

============= =============================================================== ======= ============================================================================================================ 
Name          Type                                                            Default Description                                                                                                  
============= =============================================================== ======= ============================================================================================================ 
maxCoarseDof  globalIndex                                                     0       limit of coarsening across all ranks (i.e. trim the grid hierarchy globally)                                 
partitionType geos_LinearSolverParameters_Multiscale_Coarsening_PartitionType graph   Partition type for generating coarse aggregates, valid options: ``graph``, ``cartesian``, ``semistructured`` 
ratio         real64_array                                                    {0}     Coarsening ratio (number of fine cells per coarse cell)                                                      
Graph         node                                                            unique  :ref:`XML_Graph`                                                                                             
Structured    node                                                            unique  :ref:`XML_Structured`                                                                                        
============= =============================================================== ======= ============================================================================================================ 


