

========================= ============ ================================ =========================================================================== 
Name                      Type         Registered On                    Description                                                                 
========================= ============ ================================ =========================================================================== 
averagePressure           real64                                        (no description available)                                                  
averageTemperature        real64                                        (no description available)                                                  
maxStableDt               real64                                        Value of the Maximum Stable Timestep for this solver.                       
maximumPressure           real64                                        (no description available)                                                  
maximumTemperature        real64                                        (no description available)                                                  
minimumPressure           real64                                        (no description available)                                                  
minimumTemperature        real64                                        (no description available)                                                  
phasePoreVolume           real64_array                                  (no description available)                                                  
thermalConductivityNames  string_array                                  Name of the thermal conductivity constitutive model to use                  
totalPoreVolume           real64                                        (no description available)                                                  
deltaFacePressure         real64_array :ref:`DATASTRUCTURE_FaceManager` Accumulated face pressure updates                                           
facePressure              real64_array :ref:`DATASTRUCTURE_FaceManager` Face pressure                                                               
gravityCoefficient        real64_array :ref:`DATASTRUCTURE_FaceManager` Gravity coefficient (dot product of gravity acceleration by gravity vector) 
mimGravityCoefficient     real64_array :ref:`DATASTRUCTURE_FaceManager` Mimetic gravity coefficient                                                 
LinearSolverParameters    node                                          :ref:`DATASTRUCTURE_LinearSolverParameters`                                 
NonlinearSolverParameters node                                          :ref:`DATASTRUCTURE_NonlinearSolverParameters`                              
========================= ============ ================================ =========================================================================== 


