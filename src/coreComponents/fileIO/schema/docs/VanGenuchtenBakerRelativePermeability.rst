

========================== ============ ======== ==================================================================================================== 
Name                       Type         Default  Description                                                                                          
========================== ============ ======== ==================================================================================================== 
gasOilRelPermExponentInv   real64_array 0.5      Rel perm power law exponent inverse for the pair (gas phase, oil phase) at residual water saturation 
gasOilRelPermMaxValue      real64_array 0        Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation              
name                       string       required A name is required for any non-unique nodes                                                          
phaseMinVolumeFraction     real64_array 0        Minimum volume fraction value for each phase                                                         
phaseNames                 string_array required List of fluid phases                                                                                 
waterOilRelPermExponentInv real64_array 0.5      Rel perm power law exponent inverse for the pair (water phase, oil phase) at residual gas saturation 
waterOilRelPermMaxValue    real64_array 0        Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation              
========================== ============ ======== ==================================================================================================== 


