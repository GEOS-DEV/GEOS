

============ ============ ======== ===================================================================== 
Name         Type         Default  Description                                                           
============ ============ ======== ===================================================================== 
dimensions   real64_array required Length and width of the bounded plane                                 
lengthVector R1Tensor     required Tangent vector defining the orthonormal basis along with the normal.  
name         string       required A name is required for any non-unique nodes                           
normal       R1Tensor     required Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)  
origin       R1Tensor     required Origin point (x,y,z) of the plane (basically, any point on the plane) 
widthVector  R1Tensor     required Tangent vector defining the orthonormal basis along with the normal.  
============ ============ ======== ===================================================================== 


