

========================= ============================================================================================================================== ================================================ ================================================================ 
Name                      Type                                                                                                                           Registered On                                    Description                                                      
========================= ============================================================================================================================== ================================================ ================================================================ 
maxStableDt               real64                                                                                                                                                                          Value of the Maximum Stable Timestep for this solver.            
meshTargets               geos_mapBase< std_pair< string, string >, std_vector< string, std_allocator< string > >, std_integral_constant< bool, true > >                                                  MeshBody/Region combinations that the solver will be applied to. 
parentEdgeIndex           integer_array                                                                                                                  :ref:`DATASTRUCTURE_embeddedSurfacesNodeManager` Index of parent edge within the mesh object it is registered on. 
LinearSolverParameters    node                                                                                                                                                                            :ref:`DATASTRUCTURE_LinearSolverParameters`                      
NonlinearSolverParameters node                                                                                                                                                                            :ref:`DATASTRUCTURE_NonlinearSolverParameters`                   
SolverStatistics          node                                                                                                                                                                            :ref:`DATASTRUCTURE_SolverStatistics`                            
========================= ============================================================================================================================== ================================================ ================================================================ 


