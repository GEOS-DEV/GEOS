

<<<<<<< HEAD
========================= ====== ===================================================== 
Name                      Type   Description                                           
========================= ====== ===================================================== 
maxStableDt               real64 Value of the Maximum Stable Timestep for this solver. 
LinearSolverParameters    node   :ref:`DATASTRUCTURE_LinearSolverParameters`           
NonlinearSolverParameters node   :ref:`DATASTRUCTURE_NonlinearSolverParameters`        
========================= ====== ===================================================== 
=======
========================= ================ ================================ ======================================================================================================================================================================================================================================================================================================================== 
Name                      Type             Registered On                    Description                                                                                                                                                                                                                                                                                                              
========================= ================ ================================ ======================================================================================================================================================================================================================================================================================================================== 
discretization            string                                            Name of discretization object (defined in the :ref:`NumericalMethodsManager`) to use for this solver. For instance, if this is a Finite Element Solver, the name of a :ref:`FiniteElement` should be specified. If this is a Finite Volume Method, the name of a :ref:`FiniteVolume` discretization should be specified. 
maxStableDt               real64                                            Value of the Maximum Stable Timestep for this solver.                                                                                                                                                                                                                                                                    
childIndex                localIndex_array :ref:`DATASTRUCTURE_edgeManager` Index of child within the mesh object it is registered on.                                                                                                                                                                                                                                                               
parentIndex               localIndex_array :ref:`DATASTRUCTURE_edgeManager` Index of parent within the mesh object it is registered on.                                                                                                                                                                                                                                                              
LinearSolverParameters    node                                              :ref:`DATASTRUCTURE_LinearSolverParameters`                                                                                                                                                                                                                                                                              
NonlinearSolverParameters node                                              :ref:`DATASTRUCTURE_NonlinearSolverParameters`                                                                                                                                                                                                                                                                           
========================= ================ ================================ ======================================================================================================================================================================================================================================================================================================================== 
>>>>>>> origin/develop


