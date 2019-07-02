

==================== ======= ======== ====================================================================== 
Name                 Type    Default  Description                                                            
==================== ======= ======== ====================================================================== 
beginTime            real64  0        Start time of this event                                               
endTime              real64  1e+100   End time of this event                                                 
forceDt              real64  -1       Forced timestep for this event                                         
maxEventDt           real64  -1       Forced timestep for this event                                         
maxRuntime           real64  required max runtime                                                            
name                 string  required A name is required for any non-unique nodes                            
target               string           event target                                                           
targetExactStartStop integer 1        allows timesteps to be truncated to match the start/stop times exactly 
verbosity            integer 0        Verbosity level                                                        
HaltEvent            node             :ref:`XML_HaltEvent`                                                   
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                               
SoloEvent            node             :ref:`XML_SoloEvent`                                                   
==================== ======= ======== ====================================================================== 


