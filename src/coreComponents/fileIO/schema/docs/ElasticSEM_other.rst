

========================== ================== ================================ ======================================================================= 
Name                       Type               Registered On                    Description                                                             
========================== ================== ================================ ======================================================================= 
displacementNp1AtReceivers real64_array2d                                      Displacement value at each receiver for each timestep                   
maxStableDt                real64                                              Value of the Maximum Stable Timestep for this solver.                   
receiverIsLocal            localIndex_array                                    Flag that indicates whether the receiver is local to this MPI rank      
receiverNodeIds            localIndex_array2d                                  Indices of the nodes (in the right order) for each receiver point       
sourceConstants            real64_array2d                                      Constant part of the receiver for the nodes listed in m_receiverNodeIds 
sourceIsLocal              localIndex_array                                    Flag that indicates whether the source is local to this MPI rank        
sourceNodeIds              localIndex_array2d                                  Indices of the nodes (in the right order) for each source point         
dampingVector              real64_array       :ref:`DATASTRUCTURE_nodeManager` Diagonal of the Damping Matrix.                                         
displacementx_n            real64_array       :ref:`DATASTRUCTURE_nodeManager` x-component of displacement at time n.                                  
displacementx_nm1          real64_array       :ref:`DATASTRUCTURE_nodeManager` x-component of displacement at time n-1.                                
displacementx_np1          real64_array       :ref:`DATASTRUCTURE_nodeManager` x-component of displacement at time n+1.                                
displacementy_n            real64_array       :ref:`DATASTRUCTURE_nodeManager` y-component of displacement at time n.                                  
displacementy_nm1          real64_array       :ref:`DATASTRUCTURE_nodeManager` y-component of displacement at time n-1.                                
displacementy_np1          real64_array       :ref:`DATASTRUCTURE_nodeManager` y-component of displacement at time n+1.                                
displacementz_n            real64_array       :ref:`DATASTRUCTURE_nodeManager` z-component of displacement at time n.                                  
displacementz_nm1          real64_array       :ref:`DATASTRUCTURE_nodeManager` z-component of displacement at time n-1.                                
displacementz_np1          real64_array       :ref:`DATASTRUCTURE_nodeManager` z-component of displacement at time n+1.                                
freeSurfaceFaceIndicator   localIndex_array   :ref:`DATASTRUCTURE_FaceManager` Free surface indicator, 1 if a face is on free surface 0 otherwise.     
freeSurfaceNodeIndicator   localIndex_array   :ref:`DATASTRUCTURE_nodeManager` Free surface indicator, 1 if a node is on free surface 0 otherwise.     
massVector                 real64_array       :ref:`DATASTRUCTURE_nodeManager` Diagonal of the Mass Matrix.                                            
rhs                        real64_array       :ref:`DATASTRUCTURE_nodeManager` RHS                                                                     
stiffnessVector_x          real64_array       :ref:`DATASTRUCTURE_nodeManager` x-component of stiffness vector.                                        
stiffnessVector_y          real64_array       :ref:`DATASTRUCTURE_nodeManager` y-component of stiffness vector.                                        
stiffnessVector_z          real64_array       :ref:`DATASTRUCTURE_nodeManager` z-component of stiffness vector.                                        
LinearSolverParameters     node                                                :ref:`DATASTRUCTURE_LinearSolverParameters`                             
NonlinearSolverParameters  node                                                :ref:`DATASTRUCTURE_NonlinearSolverParameters`                          
========================== ================== ================================ ======================================================================= 


