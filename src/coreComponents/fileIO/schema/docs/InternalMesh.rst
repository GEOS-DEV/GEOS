

=============== ============= ======== ================================================================================== 
Name            Type          Default  Description                                                                        
=============== ============= ======== ================================================================================== 
xCoords         real64_array  required x-coordinates of each mesh block vertex                                            
yCoords         real64_array  required y-coordinates of each mesh block vertex                                            
zCoords         real64_array  required z-coordinates of each mesh block vertex                                            
nx              integer_array required number of elements in the x-direction within each mesh block                       
ny              integer_array required number of elements in the y-direction within each mesh block                       
nz              integer_array required number of elements in the z-direction within each mesh block                       
xbias           real64_array  1        (no description available)                                                         
ybias           real64_array  1        (no description available)                                                         
zbias           real64_array  1        (no description available)                                                         
cellBlockNames  string_array  required names of each mesh block                                                           
elementTypes    string_array  required element types of each mesh block                                                   
trianglePattern integer       0        pattern by which to decompose the hex mesh into prisms (more explanation required) 
name            string        required A name is required for any non-unique nodes                                        
=============== ============= ======== ================================================================================== 


