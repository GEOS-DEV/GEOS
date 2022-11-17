

========================= =============================================================================================================================================================== ================================================================ 
Name                      Type                                                                                                                                                            Description                                                      
========================= =============================================================================================================================================================== ================================================================ 
failCriterion             integer                                                                                                                                                         (no description available)                                       
maxStableDt               real64                                                                                                                                                          Value of the Maximum Stable Timestep for this solver.            
meshTargets               geosx_mapBase< std_pair< string, string >, LvArray_Array< string, 1, camp_int_seq< long, 0l >, int, LvArray_ChaiBuffer >, std_integral_constant< bool, true > > MeshBody/Region combinations that the solver will be applied to. 
tipEdges                  LvArray_SortedArray< int, int, LvArray_ChaiBuffer >                                                                                                             Set containing all the tip edges                                 
tipFaces                  LvArray_SortedArray< int, int, LvArray_ChaiBuffer >                                                                                                             Set containing all the tip faces                                 
tipNodes                  LvArray_SortedArray< int, int, LvArray_ChaiBuffer >                                                                                                             Set containing all the nodes at the fracture tip                 
trailingFaces             LvArray_SortedArray< int, int, LvArray_ChaiBuffer >                                                                                                             Set containing all the trailing faces                            
LinearSolverParameters    node                                                                                                                                                            :ref:`DATASTRUCTURE_LinearSolverParameters`                      
NonlinearSolverParameters node                                                                                                                                                            :ref:`DATASTRUCTURE_NonlinearSolverParameters`                   
SolverStatistics          node                                                                                                                                                            :ref:`DATASTRUCTURE_SolverStatistics`                            
========================= =============================================================================================================================================================== ================================================================ 


