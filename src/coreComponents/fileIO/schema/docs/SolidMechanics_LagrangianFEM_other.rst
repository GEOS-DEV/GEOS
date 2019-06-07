

========================= ===================== ================================ ================================================================================================================================================================ 
Name                      Type                  Registered On                    Description                                                                                                                                                      
========================= ===================== ================================ ================================================================================================================================================================ 
gravityVector             R1Tensor                                               (no description available)                                                                                                                                       
maxStableDt               real64                                                 Value of the Maximum Stable Timestep for this solver.                                                                                                            
timeIntegrationOptionEnum timeIntegrationOption                                  Time integration enum class value.                                                                                                                               
Acceleration              r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the current acceleration on the nodes. This array also is used to hold the summation of nodal forces resulting from the governing equations. 
IncrementalDisplacement   r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the incremental displacements for the current time step on the nodes.                                                                        
Mass                      real64_array          :ref:`DATASTRUCTURE_nodeManager` An array that holds the mass on the nodes.                                                                                                                       
TotalDisplacement         r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the total displacements on the nodes.                                                                                                        
Velocity                  r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the current velocity on the nodes.                                                                                                           
uhatTilde                 r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the incremental displacement predictors on the nodes.                                                                                        
velocityTilde             r1_array              :ref:`DATASTRUCTURE_nodeManager` An array that holds the velocity predictors on the nodes.                                                                                                        
SystemSolverParameters    node                                                   :ref:`DATASTRUCTURE_SystemSolverParameters`                                                                                                                      
========================= ===================== ================================ ================================================================================================================================================================ 


