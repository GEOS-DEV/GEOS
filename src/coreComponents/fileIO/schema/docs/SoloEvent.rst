

==================== ======= ======== ====================================================================== 
Name                 Type    Default  Description                                                            
==================== ======= ======== ====================================================================== 
allowSubstep         integer 0        allows event sub-stepping                                              
allowSuperstep       integer 0        allows event super-stepping (dt_super=dt+t-t_last)                     
beginTime            real64  0        Start time of this event                                               
endTime              real64  1e+100   End time of this event                                                 
forceDt              real64  -1       Forced timestep for this event                                         
name                 string  required A name is required for any non-unique nodes                            
substepFactor        integer 1        integer substep factor (dt_sub=dt/f)                                   
target               string  required event target                                                           
targetCycle          integer -1       event cycle                                                            
targetExactStartStop integer 0        allows timesteps to be truncated to match the start/stop times exactly 
targetExactTimestep  integer -1       allows timesteps to be truncated to match time frequency perfectly     
targetTime           real64  -1       Event time                                                             
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
==================== ======= ======== ====================================================================== 


