

========================= ================================= ================================ =========================================================== 
Name                      Type                              Registered On                    Description                                                 
========================= ================================= ================================ =========================================================== 
failCriterion             integer                                                            (no description available)                                  
maxStableDt               real64                                                             Value of the Maximum Stable Timestep for this solver.       
tipEdges                  LvArray_SortedArray< long, long >                                  Set containing all the tip edges                            
tipFaces                  LvArray_SortedArray< long, long >                                  Set containing all the tip faces                            
tipNodes                  LvArray_SortedArray< long, long >                                  Set containing all the nodes at the fracture tip            
trailingFaces             LvArray_SortedArray< long, long >                                  Set containing all the trailing faces                       
SIFNode                   real64_array                      :ref:`DATASTRUCTURE_nodeManager` SIF on the node.                                            
SIF_I                     real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_I of the edge.                                          
SIF_II                    real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_II of the edge.                                         
SIF_III                   real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_III of the edge.                                        
childIndex                localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Index of child within the  mesh object it is registered on. 
degreeFromCrack           integer_array                     :ref:`DATASTRUCTURE_nodeManager` Connectivity distance from crack.                           
degreeFromCrackTip        integer_array                     :ref:`DATASTRUCTURE_nodeManager` Degree of connectivity separation from crack tip.           
parentIndex               localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Index of parent within the mesh object it is registered on. 
ruptureTime               real64_array                      :ref:`DATASTRUCTURE_nodeManager` Time that the node was ruptured.                            
LinearSolverParameters    node                                                               :ref:`DATASTRUCTURE_LinearSolverParameters`                 
NonlinearSolverParameters node                                                               :ref:`DATASTRUCTURE_NonlinearSolverParameters`              
========================= ================================= ================================ =========================================================== 


