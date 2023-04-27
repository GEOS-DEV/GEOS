

============ ======== ======== =================================================================================================================== 
Name         Type     Default  Description                                                                                                         
============ ======== ======== =================================================================================================================== 
center       R1Tensor required (x,y,z) coordinates of the center of the RTheta                                                                     
lengthVector R1Tensor required Tangent vector defining the orthonormal basis along with the normal.                                                
name         string   required A name is required for any non-unique nodes                                                                         
normal       R1Tensor required Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)                                                
tolerance    real64   1e-05    Tolerance to determine if a point sits on the RTheta or not. It is relative to the maximum dimension of the RTheta. 
widthVector  R1Tensor required Tangent vector defining the orthonormal basis along with the normal.                                                
============ ======== ======== =================================================================================================================== 


