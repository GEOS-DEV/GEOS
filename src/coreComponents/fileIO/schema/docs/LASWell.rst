

====================== ======= ======== ====================================================================================== 
Name                   Type    Default  Description                                                                            
====================== ======= ======== ====================================================================================== 
fileName               path    required Path to the las file                                                                   
geometryLogIndexInFile integer -1       Position of the log to take if there are several log sections defined in the LAS file  
meshName               string  required Name of the reservoir mesh associated with this well                                   
name                   string  required A name is required for any non-unique nodes                                            
numElementsPerSegment  integer required Number of well elements per polyline segment                                           
radius                 real64  required Radius of the well                                                                     
wellControlsName       string  required Name of the set of constraints associated with this well                               
wellRegionName         string  required Name of the well element region                                                        
Perforation            node             :ref:`XML_Perforation`                                                                 
====================== ======= ======== ====================================================================================== 


