

========================= ============ ================================ =========================================================================== 
Name                      Type         Registered On                    Description                                                                 
========================= ============ ================================ =========================================================================== 
maxStableDt               real64                                        Value of the Maximum Stable Timestep for this solver.                       
deltaFacePressure         real64_array :ref:`DATASTRUCTURE_FaceManager` Accumulated face pressure updates                                           
facePressure              real64_array :ref:`DATASTRUCTURE_FaceManager` Face pressure                                                               
gravityCoefficient        real64_array :ref:`DATASTRUCTURE_FaceManager` Gravity coefficient (dot product of gravity acceleration by gravity vector) 
LinearSolverParameters    node                                          :ref:`DATASTRUCTURE_LinearSolverParameters`                                 
NonlinearSolverParameters node                                          :ref:`DATASTRUCTURE_NonlinearSolverParameters`                              
========================= ============ ================================ =========================================================================== 


