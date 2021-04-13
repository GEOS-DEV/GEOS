

========================= ================== ================================ ======================================================================= 
Name                      Type               Registered On                    Description                                                             
========================= ================== ================================ ======================================================================= 
maxStableDt               real64                                              Value of the Maximum Stable Timestep for this solver.                   
pressureNp1AtReceivers    real64_array                                        Pressure value at each receiver for each timestep                       
receiverIsLocal           localIndex_array                                    Flag that indicates whether the receiver is local to this MPI rank      
receiverNodeIds           localIndex_array2d                                  Indices of the nodes (in the right order) for each receiver point       
sourceConstants           real64_array2d                                      Constant part of the receiver for the nodes listed in m_receiverNodeIds 
sourceIsLocal             localIndex_array                                    Flag that indicates whether the source is local to this MPI rank        
sourceNodeIds             localIndex_array2d                                  Indices of the nodes (in the right order) for each source point         
dampingVector             real64_array       :ref:`DATASTRUCTURE_nodeManager` Diagonal of the Damping Matrix.                                         
freeSurfaceFaceIndicator  localIndex_array   :ref:`DATASTRUCTURE_FaceManager` Free surface indicator, 1 if a face is on free surface 0 otherwise.     
freeSurfaceNodeIndicator  localIndex_array   :ref:`DATASTRUCTURE_nodeManager` Free surface indicator, 1 if a node is on free surface 0 otherwise.     
massVector                real64_array       :ref:`DATASTRUCTURE_nodeManager` Diagonal of the Mass Matrix.                                            
pressure_n                real64_array       :ref:`DATASTRUCTURE_nodeManager` Scalar pressure at time n.                                              
pressure_nm1              real64_array       :ref:`DATASTRUCTURE_nodeManager` Scalar pressure at time n-1.                                            
pressure_np1              real64_array       :ref:`DATASTRUCTURE_nodeManager` Scalar pressure at time n+1.                                            
rhs                       real64_array       :ref:`DATASTRUCTURE_nodeManager` RHS                                                                     
stiffnessVector           real64_array       :ref:`DATASTRUCTURE_nodeManager` Stiffness vector contains R_h*Pressure_n.                               
LinearSolverParameters    node                                                :ref:`DATASTRUCTURE_LinearSolverParameters`                             
NonlinearSolverParameters node                                                :ref:`DATASTRUCTURE_NonlinearSolverParameters`                          
========================= ================== ================================ ======================================================================= 


