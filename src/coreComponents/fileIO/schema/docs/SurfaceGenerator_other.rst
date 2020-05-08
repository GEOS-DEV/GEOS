

========================= ================================= ================================ ===================================================== 
Name                      Type                              Registered On                    Description                                           
========================= ================================= ================================ ===================================================== 
failCriterion             integer                                                            (no description available)                            
maxStableDt               real64                                                             Value of the Maximum Stable Timestep for this solver. 
tipEdges                  LvArray_SortedArray< long, long >                                  Set containing all the tip edges                      
tipFaces                  LvArray_SortedArray< long, long >                                  Set containing all the tip faces                      
tipNodes                  LvArray_SortedArray< long, long >                                  Set containing all the nodes at the fracture tip      
trailingFaces             LvArray_SortedArray< long, long >                                  Set containing all the trailing faces                 
childIndex                localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Child index of the edge.                              
parentIndex               localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Parent index of the edge.                             
NonlinearSolverParameters node                                                               :ref:`DATASTRUCTURE_NonlinearSolverParameters`        
SystemSolverParameters    node                                                               :ref:`DATASTRUCTURE_SystemSolverParameters`           
========================= ================================= ================================ ===================================================== 


