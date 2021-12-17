

================== ============ ======== ============================================= 
Name               Type         Default  Description                                   
================== ============ ======== ============================================= 
baseline           path         none     Baseline file                                 
feedComposition    real64_array required Feed composition array [mol fraction]         
fluid              string       required Fluid to test                                 
logLevel           integer      0        Log level                                     
name               string       required A name is required for any non-unique nodes   
output             string       none     Output file                                   
pressureControl    string       required Function controlling pressure time history    
steps              integer      required Number of load steps to take                  
temperatureControl string       required Function controlling temperature time history 
================== ============ ======== ============================================= 


