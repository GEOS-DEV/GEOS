

======================== ============ ======== =========================================== 
Name                     Type         Default  Description                                 
======================== ============ ======== =========================================== 
componentNames           string_array {}       List of fluid component names               
compressibility          real64       0        Fluid compressibility                       
defaultCompressibility   real64_array {0}      Default value for compressibility.          
defaultDensity           real64_array {0}      Default value for density.                  
defaultViscosity         real64_array {0}      Default value for viscosity.                
flowBehaviorIndex        real64_array {0}      Flow behavior index                         
flowConsistencyIndex     real64_array {0}      Flow consistency index                      
maxProppantConcentration real64       0.6      Maximum proppant concentration              
name                     string       required A name is required for any non-unique nodes 
referenceDensity         real64       1000     Reference fluid density                     
referencePressure        real64       100000   Reference pressure                          
referenceProppantDensity real64       1400     Reference proppant density                  
referenceViscosity       real64       0.001    Reference fluid viscosity                   
======================== ============ ======== =========================================== 


