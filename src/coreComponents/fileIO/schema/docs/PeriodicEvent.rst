

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
timeFrequency        real64  -1       event frequency (time)                                                 
cycleFrequency       integer 1        event frequency (cycle, Default)                                       
targetExactTimestep  integer -1       allows timesteps to be truncated to match time frequency perfectly     
function             string           Name of the symbolic math function                                     
object               string           Path of the function input object (directory format)                   
set                  string           Setname of the input object (if empty, default to everything)          
stat                 integer 0        Selection of the min/avg/max for functions that target vectors         
threshold            real64  1e+10    event threshold                                                        
name                 string  required A name is required for any non-unique nodes                            
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
==================== ======= ======== ====================================================================== 


