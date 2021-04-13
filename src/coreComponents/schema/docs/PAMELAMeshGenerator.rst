

================= ============ ======== ======================================================== 
Name              Type         Default  Description                                              
================= ============ ======== ======================================================== 
fieldNamesInGEOSX string_array {}       Name of the fields within GEOSX                          
fieldsToImport    string_array {}       Fields to be imported from the external mesh file        
file              path         required path to the mesh file                                    
name              string       required A name is required for any non-unique nodes              
reverseZ          integer      0        0 : Z coordinate is upward, 1 : Z coordinate is downward 
scale             real64       1        Scale the coordinates of the vertices                    
================= ============ ======== ======================================================== 


