

=============== ============================================================== ======= ========================================================================================= 
Name            Type                                                           Default Description                                                                               
=============== ============================================================== ======= ========================================================================================= 
matrixWeights   integer                                                        0       If >0, specifies matrix weight multiplier when building graph edge weights                
method          geos_LinearSolverParameters_Multiscale_Coarsening_Graph_Method metis   Graph partitioning method, valid options: ``metis``, ``scotch``                           
minCommonNodes  integer                                                        3       Minimum number of nodes shared between two cells when constructing the connectivity graph 
preserveRegions integer                                                        0       Attempt to keep cells from the same region in one aggregate                               
Metis           node                                                           unique  :ref:`XML_Metis`                                                                          
=============== ============================================================== ======= ========================================================================================= 


