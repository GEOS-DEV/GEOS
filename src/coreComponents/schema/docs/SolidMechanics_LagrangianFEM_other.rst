

========================= ============== ================================ ================================================================================================================================================================ 
Name                      Type           Registered On                    Description                                                                                                                                                      
========================= ============== ================================ ================================================================================================================================================================ 
maxForce                  real64                                          The maximum force contribution in the problem domain.                                                                                                            
maxStableDt               real64                                          Value of the Maximum Stable Timestep for this solver.                                                                                                            
Acceleration              real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the current acceleration on the nodes. This array also is used to hold the summation of nodal forces resulting from the governing equations. 
IncrementalDisplacement   real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the incremental displacements for the current time step on the nodes.                                                                        
Mass                      real64_array   :ref:`DATASTRUCTURE_nodeManager` An array that holds the mass on the nodes.                                                                                                                       
TotalDisplacement         real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the total displacements on the nodes.                                                                                                        
Velocity                  real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the current velocity on the nodes.                                                                                                           
contactForce              real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the contact force.                                                                                                                           
externalForce             real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the external forces on the nodes. This includes any boundary conditions as well as coupling forces such as hydraulic forces.                 
uhatTilde                 real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the incremental displacement predictors on the nodes.                                                                                        
velocityTilde             real64_array2d :ref:`DATASTRUCTURE_nodeManager` An array that holds the velocity predictors on the nodes.                                                                                                        
LinearSolverParameters    node                                            :ref:`DATASTRUCTURE_LinearSolverParameters`                                                                                                                      
NonlinearSolverParameters node                                            :ref:`DATASTRUCTURE_NonlinearSolverParameters`                                                                                                                   
========================= ============== ================================ ================================================================================================================================================================ 


