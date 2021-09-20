

========================= ============== ================================ ================================================================================================================================================================ 
Name                      Type           Registered On                    Description                                                                                                                                                      
========================= ============== ================================ ================================================================================================================================================================ 
maxStableDt               real64                                          Value of the Maximum Stable Timestep for this solver.                                                                                                            
Acceleration              real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the current acceleration on the nodes. This array also is used to hold the summation of nodal forces resulting from the governing equations. 
Mass                      real64_array   :ref:`DATASTRUCTURE_nodeManager` An array that holds the mass on the nodes.                                                                                                                       
Velocity                  real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the current velocity on the nodes.                                                                                                           
LinearSolverParameters    node                                            :ref:`DATASTRUCTURE_LinearSolverParameters`                                                                                                                      
NonlinearSolverParameters node                                            :ref:`DATASTRUCTURE_NonlinearSolverParameters`                                                                                                                   
========================= ============== ================================ ================================================================================================================================================================ 


