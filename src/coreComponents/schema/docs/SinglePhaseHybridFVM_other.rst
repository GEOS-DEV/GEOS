

========================= ==================================================================================================================================================== ================================ ================================================================ 
Name                      Type                                                                                                                                                 Registered On                    Description                                                      
========================= ==================================================================================================================================================== ================================ ================================================================ 
maxStableDt               real64                                                                                                                                                                                Value of the Maximum Stable Timestep for this solver.            
meshTargets               geosx_mapBase< std_string, LvArray_Array< std_string, 1, camp_int_seq< long, 0l >, long, LvArray_ChaiBuffer >, std_integral_constant< bool, true > >                                  MeshBody/Region combinations that the solver will be applied to. 
deltaFacePressure         real64_array                                                                                                                                         :ref:`DATASTRUCTURE_FaceManager` Accumulated face pressure updates                                
LinearSolverParameters    node                                                                                                                                                                                  :ref:`DATASTRUCTURE_LinearSolverParameters`                      
NonlinearSolverParameters node                                                                                                                                                                                  :ref:`DATASTRUCTURE_NonlinearSolverParameters`                   
========================= ==================================================================================================================================================== ================================ ================================================================ 


