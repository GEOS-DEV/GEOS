

==================== ==================================== ======== ============================================================= 
Name                 Type                                 Default  Description                                                   
==================== ==================================== ======== ============================================================= 
componentMolarWeight real64_array                         required Component molar weights                                       
componentNames       string_array                         {}       List of component names                                       
fluidType            constitutive_BlackOilFluid_FluidType required | Type of black-oil fluid. Valid options:                       
                                                                   | * DeadOil                                                     
                                                                   | * LiveOil                                                     
name                 string                               required A name is required for any non-unique nodes                   
phaseNames           string_array                         required List of fluid phases                                          
surfaceDensities     real64_array                         required List of surface densities for each phase                      
tableFiles           Path_array                           required List of filenames with input PVT tables                       
==================== ==================================== ======== ============================================================= 


