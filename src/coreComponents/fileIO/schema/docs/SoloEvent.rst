

==================== ======= ======== ===================================================================================================================== 
Name                 Type    Default  Description                                                                                                           
==================== ======= ======== ===================================================================================================================== 
beginTime            real64  0        Start time of this event.                                                                                             
endTime              real64  1e+100   End time of this event.                                                                                               
forceDt              real64  -1       While active, this event will request this timestep value (ignoring any children/targets requests).                   
maxEventDt           real64  -1       While active, this event will request a timestep <= this value (depending upon any child/target requests).            
name                 string  required A name is required for any non-unique nodes                                                                           
target               string           Name of the object to be executed when the event criteria are met.                                                    
targetCycle          integer -1       Target event cycle.                                                                                                   
targetExactStartStop integer 1        If this option is set, the event will reduce its timestep requests to match any specified beginTime/endTimes exactly. 
targetExactTimestep  integer 1        If this option is set, the event will reduce its timestep requests to match the specified execution time exactly.     
targetTime           real64  -1       Target event time.                                                                                                    
verbosity            integer 0        Verbosity level                                                                                                       
HaltEvent            node             :ref:`XML_HaltEvent`                                                                                                  
PeriodicEvent        node             :ref:`XML_PeriodicEvent`                                                                                              
SoloEvent            node             :ref:`XML_SoloEvent`                                                                                                  
==================== ======= ======== ===================================================================================================================== 


