

========================= ============================================================================================================================== ======================================================================= 
Name                      Type                                                                                                                           Description                                                             
========================= ============================================================================================================================== ======================================================================= 
indexSeismoTrace          integer                                                                                                                        Count for output pressure at receivers                                  
maxStableDt               real64                                                                                                                         Value of the Maximum Stable Timestep for this solver.                   
meshTargets               geos_mapBase< std_pair< string, string >, std_vector< string, std_allocator< string > >, std_integral_constant< bool, true > > MeshBody/Region combinations that the solver will be applied to.        
pressureNp1AtReceivers    real32_array2d                                                                                                                 Pressure value at each receiver for each timestep                       
rcvElem                   integer_array                                                                                                                  Element containing the receivers                                        
receiverConstants         real64_array2d                                                                                                                 Constant part of the receiver for the nodes listed in m_receiverNodeIds 
receiverIsLocal           integer_array                                                                                                                  Flag that indicates whether the receiver is local to this MPI rank      
receiverNodeIds           integer_array2d                                                                                                                Indices of the nodes (in the right order) for each receiver point       
receiverRegion            integer_array                                                                                                                  Region containing the receivers                                         
sourceConstants           real64_array2d                                                                                                                 Constant part of the source for the nodes listed in m_sourceNodeIds     
sourceIsAccessible        integer_array                                                                                                                  Flag that indicates whether the source is local to this MPI rank        
sourceNodeIds             integer_array2d                                                                                                                Indices of the nodes (in the right order) for each source point         
sourceValue               real32_array2d                                                                                                                 Source Value of the sources                                             
useDAS                    integer                                                                                                                        Flag to indicate if DAS type of data will be modeled                    
usePML                    integer                                                                                                                        Flag to apply PML                                                       
LinearSolverParameters    node                                                                                                                           :ref:`DATASTRUCTURE_LinearSolverParameters`                             
NonlinearSolverParameters node                                                                                                                           :ref:`DATASTRUCTURE_NonlinearSolverParameters`                          
SolverStatistics          node                                                                                                                           :ref:`DATASTRUCTURE_SolverStatistics`                                   
========================= ============================================================================================================================== ======================================================================= 


