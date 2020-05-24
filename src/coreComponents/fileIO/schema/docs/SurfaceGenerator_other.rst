

========================= ================================= ======================================================================================================== 
Name                      Type                              Description                                                                                              
========================= ================================= ======================================================================================================== 
failCriterion             integer                           (no description available)                                                                               
maxStableDt               real64                            Value of the Maximum Stable Timestep for this solver.                                                    
nodesWithAssignedDisp     LvArray_SortedArray< long, long > Set containing all the nodes with displacement assigned due to volume in the partially fractured element 
tipEdges                  LvArray_SortedArray< long, long > Set containing all the tip edges                                                                         
tipFaces                  LvArray_SortedArray< long, long > Set containing all the tip faces                                                                         
tipNodes                  LvArray_SortedArray< long, long > Set containing all the nodes at the fracture tip                                                         
trailingFaces             LvArray_SortedArray< long, long > Set containing all the trailing faces                                                                    
LinearSolverParameters    node                              :ref:`DATASTRUCTURE_LinearSolverParameters`                                                              
NonlinearSolverParameters node                              :ref:`DATASTRUCTURE_NonlinearSolverParameters`                                                           
========================= ================================= ======================================================================================================== 


