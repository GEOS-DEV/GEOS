

=============== ======= ================================== ============= 
Name            Type    Description                        Registered By 
=============== ======= ================================== ============= 
time            real64  Current simulation time.                         
dt              real64  Current simulation timestep.                     
cycle           integer Current simulation cycle number.                 
currentSubEvent integer index of the current subevent.                   
currentMaxDt    real64  Maximum dt request for event loop.               
HaltEvent       node    :ref:`DATASTRUCTURE_HaltEvent`                   
SoloEvent       node    :ref:`DATASTRUCTURE_SoloEvent`                   
PeriodicEvent   node    :ref:`DATASTRUCTURE_PeriodicEvent`               
=============== ======= ================================== ============= 


