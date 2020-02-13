

===================== =================== ======== ======================================================== 
Name                  Type                Default  Description                                              
===================== =================== ======== ======================================================== 
crossSectionArea      real64              required cross section area of the well                           
meshName              string              required name of the reservoir mesh associated with this well     
name                  string              required A name is required for any non-unique nodes              
numElementsPerSegment integer             required number of well elements per polyline segment             
polylineNodeCoords    real64_array2d      required physical coordinates of the well polyline nodes          
polylineSegmentConn   globalIndex_array2d required connectivity of the polyline segments                    
wellControlsName      string              required name of the set of constraints associated with this well 
wellRegionName        string              required name of the well element region                          
Perforation           node                         :ref:`XML_Perforation`                                   
===================== =================== ======== ======================================================== 


