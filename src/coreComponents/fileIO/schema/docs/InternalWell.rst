

========================== =================== ======== ==================================================== 
Name                       Type                Default  Description                                          
========================== =================== ======== ==================================================== 
meshName                   string              required name of the reservoir mesh associated with this well 
name                       string              required A name is required for any non-unique nodes          
nodeCoords                 real64_array2d      required physical coordinates of the well polyline            
numberOfElementsPerSegment integer             required number of well elements per polyline segment         
segmentConn                globalIndex_array2d required connectivity of the polyline segment                 
wellRegionName             string              required name of the well element region                      
========================== =================== ======== ==================================================== 


