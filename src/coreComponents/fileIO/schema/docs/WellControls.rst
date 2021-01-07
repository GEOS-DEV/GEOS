

==================== ========================== ======== =================================================================================================================================================== 
Name                 Type                       Default  Description                                                                                                                                         
==================== ========================== ======== =================================================================================================================================================== 
control              geosx_WellControls_Control required | Well control. Valid options:                                                                                                                        
                                                         | * BHP                                                                                                                                               
                                                         | * oilVolRate                                                                                                                                        
                                                         | * totalVolRate                                                                                                                                      
injectionStream      real64_array               {-1}     Global component densities for the injection stream                                                                                                 
name                 string                     required A name is required for any non-unique nodes                                                                                                         
referenceElevation   real64                     required Reference elevation where BHP control is enforced                                                                                                   
targetBHP            real64                     required Target bottom-hole pressure                                                                                                                         
targetOilRate        real64                     1e+09    Target oil rate                                                                                                                                     
targetRate           real64                     required Target rate                                                                                                                                         
type                 geosx_WellControls_Type    required | Well type. Valid options:                                                                                                                           
                                                         | * producer                                                                                                                                          
                                                         | * injector                                                                                                                                          
useSurfaceConditions integer                    0        | Flag to specify whether rates are checked at surface or reservoir conditions.                                                                       
                                                         | Equal to 1 for surface conditions, and to 0 for reservoir conditions                                                                                
==================== ========================== ======== =================================================================================================================================================== 


