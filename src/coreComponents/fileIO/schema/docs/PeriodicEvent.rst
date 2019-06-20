

==================== ======= ======== ====================================================================== 
Name                 Type    Default  Description                                                            
==================== ======= ======== ====================================================================== 
allowSubstep         integer 0        allows event sub-stepping                                              
allowSuperstep       integer 0        allows event super-stepping (dt_super=dt+t-t_last)                     
beginTime            real64  0        Start time of this event                                               
cycleFrequency       integer 1        event frequency (cycle, Default)                                       
endTime              real64  1e+100   End time of this event                                                 
forceDt              real64  -1       Forced timestep for this event                                         
function             string           Name of the symbolic math function                                     
name                 string  required A name is required for any non-unique nodes                            
object               string           Path of the function input object (directory format)                   
set                  string           Setname of the input object (if empty, default to everything)          
stat                 integer 0        Selection of the min/avg/max for functions that target vectors         
substepFactor        integer 1        integer substep factor (dt_sub=dt/f)                                   
target               string  required event target                                                           
targetExactStartStop integer 0        allows timesteps to be truncated to match the start/stop times exactly 
targetExactTimestep  integer -1       allows timesteps to be truncated to match time frequency perfectly     
threshold            real64  1e+10    event threshold                                                        
timeFrequency        real64  -1       event frequency (time)                                                 
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
==================== ======= ======== ====================================================================== 


