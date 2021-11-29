

=============== =============================================== =========== =================================================================================== 
Name            Type                                            Default     Description                                                                         
=============== =============================================== =========== =================================================================================== 
defaultAperture real64                                          required    The default aperture of newly formed surface elements.                              
materialList    string_array                                    required    List of materials present in this region                                            
name            string                                          required    A name is required for any non-unique nodes                                         
subRegionType   geosx_SurfaceElementRegion_SurfaceSubRegionType faceElement | Type of surface element subregion. Valid options:                                   
                                                                            | * faceElement                                                                       
                                                                            | * embeddedElement                                                                   
=============== =============================================== =========== =================================================================================== 


