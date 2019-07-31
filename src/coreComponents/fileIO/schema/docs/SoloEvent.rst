

==================== ======= ======== ====================================================================== 
Name                 Type    Default  Description                                                            
==================== ======= ======== ====================================================================== 
beginTime            real64  0        Start time of this event                                               
endTime              real64  1e+100   End time of this event                                                 
forceDt              real64  -1       Forced timestep for this event                                         
maxEventDt           real64  -1       Forced timestep for this event                                         
name                 string  required A name is required for any non-unique nodes                            
target               string           event target                                                           
targetCycle          integer -1       Event cycle                                                            
targetExactStartStop integer 1        allows timesteps to be truncated to match the start/stop times exactly 
targetExactTimestep  integer 1        Allows timesteps to be truncated to match time frequency perfectly     
targetTime           real64  -1       Event time                                                             
verbosity            integer 0        Verbosity level                                                        
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
==================== ======= ======== ====================================================================== 


