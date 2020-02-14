

=============== ============ ======== =================================================== 
Name            Type         Default  Description                                         
=============== ============ ======== =================================================== 
control         string       required Well control (BHP/gasRate/oilRate/waterRate)        
injectionStream real64_array {-1}     Global component densities for the injection stream 
name            string       required A name is required for any non-unique nodes         
referenceDepth  real64       0        Reference depth for well bottom hole pressure       
targetBHP       real64       required Target bottom-hole pressure                         
targetRate      real64       required Target rate                                         
type            string       required Well type (producer/injector)                       
=============== ============ ======== =================================================== 


