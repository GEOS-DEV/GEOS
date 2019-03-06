

=================== ======= ================================================================== 
Name                Type    Description                                                        
=================== ======= ================================================================== 
lastTime            real64  last event occurrence (time)                                       
lastCycle           integer last event occurrence (cycle)                                      
currentSubEvent     integer index of the current subevent                                      
isTargetExecuting   integer index of the current subevent                                      
targetTime          real64  Event time                                                         
targetCycle         integer event cycle                                                        
targetExactTimestep integer allows timesteps to be truncated to match time frequency perfectly 
HaltEvent           node    :ref:`DATASTRUCTURE_HaltEvent`                                     
SoloEvent           node    :ref:`DATASTRUCTURE_SoloEvent`                                     
PeriodicEvent       node    :ref:`DATASTRUCTURE_PeriodicEvent`                                 
=================== ======= ================================================================== 


