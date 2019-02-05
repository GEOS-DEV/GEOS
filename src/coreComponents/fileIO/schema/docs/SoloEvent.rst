

==================== ======= ======== ====================================================================== 
Name                 Type    Default  Description                                                            
==================== ======= ======== ====================================================================== 
target               string  required event target                                                           
beginTime            real64  0        Start time of this event                                               
endTime              real64  1e+100   End time of this event                                                 
forceDt              real64  -1       Forced timestep for this event                                         
allowSuperstep       integer 0        allows event super-stepping (dt_super=dt+t-t_last)                     
allowSubstep         integer 0        allows event sub-stepping                                              
substepFactor        integer 1        integer substep factor (dt_sub=dt/f)                                   
targetExactStartStop integer 0        allows timesteps to be truncated to match the start/stop times exactly 
name                 string  required A name is required for any non-unique nodes                            
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
==================== ======= ======== ====================================================================== 


