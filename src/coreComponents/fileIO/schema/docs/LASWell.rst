

====================== ======= ======== ====================================================================================== 
Name                   Type    Default  Description                                                                            
====================== ======= ======== ====================================================================================== 
crossSectionArea       real64  required cross section area of the well                                                         
fileName               string  required Path to the las file                                                                   
geometryLogIndexInFile integer -1       Position of the log to take if there are several log sections defined in the LAS file  
justImportGeometry     integer 1        Just import the geometry of the well and ignore the other logs                         
meshName               string  required name of the reservoir mesh associated with this well                                   
name                   string  required A name is required for any non-unique nodes                                            
numElementsPerSegment  integer required number of well elements per polyline segment                                           
wellControlsName       string  required name of the set of constraints associated with this well                               
wellRegionName         string  required name of the well element region                                                        
====================== ======= ======== ====================================================================================== 


