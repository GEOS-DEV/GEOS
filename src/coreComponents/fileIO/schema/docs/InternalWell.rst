

===================== =================== ======== ======================================================== 
Name                  Type                Default  Description                                              
===================== =================== ======== ======================================================== 
meshName              string              required Name of the reservoir mesh associated with this well     
name                  string              required A name is required for any non-unique nodes              
numElementsPerSegment integer             required Number of well elements per polyline segment             
polylineNodeCoords    real64_array2d      required Physical coordinates of the well polyline nodes          
polylineSegmentConn   globalIndex_array2d required Connectivity of the polyline segments                    
radius                real64              required Radius of the well                                       
wellControlsName      string              required Name of the set of constraints associated with this well 
wellRegionName        string              required Name of the well element region                          
Perforation           node                         :ref:`XML_Perforation`                                   
===================== =================== ======== ======================================================== 


