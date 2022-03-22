

========================= =================================================================================================================================================== ======================================================================= 
Name                      Type                                                                                                                                                Description                                                             
========================= =================================================================================================================================================== ======================================================================= 
maxStableDt               real64                                                                                                                                              Value of the Maximum Stable Timestep for this solver.                   
meshTargets               geosx_mapBase< std_string, LvArray_Array< std_string, 1, camp_int_seq< long, 0l >, int, LvArray_ChaiBuffer >, std_integral_constant< bool, true > > MeshBody/Region combinations that the solver will be applied to.        
pressureNp1AtReceivers    real64_array                                                                                                                                        Pressure value at each receiver for each timestep                       
receiverIsLocal           integer_array                                                                                                                                       Flag that indicates whether the receiver is local to this MPI rank      
receiverNodeIds           integer_array2d                                                                                                                                     Indices of the nodes (in the right order) for each receiver point       
sourceConstants           real64_array2d                                                                                                                                      Constant part of the receiver for the nodes listed in m_receiverNodeIds 
sourceIsLocal             integer_array                                                                                                                                       Flag that indicates whether the source is local to this MPI rank        
sourceNodeIds             integer_array2d                                                                                                                                     Indices of the nodes (in the right order) for each source point         
LinearSolverParameters    node                                                                                                                                                :ref:`DATASTRUCTURE_LinearSolverParameters`                             
NonlinearSolverParameters node                                                                                                                                                :ref:`DATASTRUCTURE_NonlinearSolverParameters`                          
========================= =================================================================================================================================================== ======================================================================= 


