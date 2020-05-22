

========================= ================================= ================================ =================================================================================== 
Name                      Type                              Registered On                    Description                                                                         
========================= ================================= ================================ =================================================================================== 
failCriterion             integer                                                            (no description available)                                                          
maxStableDt               real64                                                             Value of the Maximum Stable Timestep for this solver.                               
tipEdges                  LvArray_SortedArray< long, long >                                  Set containing all the tip edges                                                    
tipFaces                  LvArray_SortedArray< long, long >                                  Set containing all the tip faces                                                    
tipNodes                  LvArray_SortedArray< long, long >                                  Set containing all the nodes at the fracture tip                                    
trailingFaces             LvArray_SortedArray< long, long >                                  Set containing all the trailing faces                                               
K_IC                      r1_array                          :ref:`DATASTRUCTURE_FaceManager` K_IC in each plane.                                                                 
SIFNode                   real64_array                      :ref:`DATASTRUCTURE_nodeManager` SIF on the node.                                                                    
SIF_I                     real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_I of the edge.                                                                  
SIF_II                    real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_II of the edge.                                                                 
SIF_III                   real64_array                      :ref:`DATASTRUCTURE_edgeManager` SIF_III of the edge.                                                                
SIFonFace                 real64_array                      :ref:`DATASTRUCTURE_FaceManager` SIF on the face.                                                                    
childIndex                localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Index of child within the  mesh object it is registered on.                         
degreeFromCrack           integer_array                     :ref:`DATASTRUCTURE_nodeManager` Connectivity distance from crack.                                                   
degreeFromCrackTip        integer_array                     :ref:`DATASTRUCTURE_nodeManager` Degree of connectivity separation from crack tip.                                   
isFaceSeparable           integer_array                     :ref:`DATASTRUCTURE_FaceManager` A flag to mark if the face is separable.                                            
parentIndex               localIndex_array                  :ref:`DATASTRUCTURE_edgeManager` Index of parent within the mesh object it is registered on.                         
primaryCandidateFace      localIndex_array                  :ref:`DATASTRUCTURE_FaceManager` SIF_III of the edge.                                                                
ruptureState              integer_array                     :ref:`DATASTRUCTURE_FaceManager` Rupture state of the face.0=not ready for rupture. 1=ready for rupture. 2=ruptured. 
ruptureTime               real64_array                      :ref:`DATASTRUCTURE_nodeManager` Time that the object was ruptured.                                                  
LinearSolverParameters    node                                                               :ref:`DATASTRUCTURE_LinearSolverParameters`                                         
NonlinearSolverParameters node                                                               :ref:`DATASTRUCTURE_NonlinearSolverParameters`                                      
========================= ================================= ================================ =================================================================================== 


