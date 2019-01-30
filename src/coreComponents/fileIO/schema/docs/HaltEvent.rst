

==================== ======= ======= ====================================================================== 
Name                 Type    Default Description                                                            
==================== ======= ======= ====================================================================== 
target               string          event target                                                           
beginTime            real64  0       Start time of this event                                               
endTime              real64  1e+100  End time of this event                                                 
forceDt              real64  -1      Forced timestep for this event                                         
allowSuperstep       integer 0       allows event super-stepping (dt_super=dt+t-t_last)                     
allowSubstep         integer 0       allows event sub-stepping                                              
substepFactor        integer 1       integer substep factor (dt_sub=dt/f)                                   
targetExactStartStop integer 0       allows timesteps to be truncated to match the start/stop times exactly 
maxRuntime           real64  0       max runtime                                                            
name                 string          A name is required for any non-unique nodes                            
HaltEvent            node            `XML_HaltEvent`_                                                       
SoloEvent            node            `XML_SoloEvent`_                                                       
PeriodicEvent        node            `XML_PeriodicEvent`_                                                   
==================== ======= ======= ====================================================================== 


