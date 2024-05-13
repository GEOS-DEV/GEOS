set( kernelPath "coreComponents/physicsSolvers/fluidFlow/kernels" )
set( fluidTypes BlackOilFluid
                CO2BrineEzrokhiFluid
                CO2BrineEzrokhiThermalFluid
                CO2BrinePhillipsFluid
                CO2BrinePhillipsThermalFluid
                CompositionalTwoPhasePengRobinsonConstantViscosity
                CompositionalTwoPhasePengRobinsonLBCViscosity
                CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity
                CompositionalTwoPhaseSoaveRedlichKwongLBCViscosity
                DeadOilFluid )

if( ENABLE_PVTPackage )
    set( fluidTypes
         ${fluidTypes}
         CompositionalMultiphaseFluidPVTPackage )
endif()

foreach( FLUID_TYPE ${fluidTypes} )
    set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/MultiFluidUpdate_${FLUID_TYPE}.cpp" )
    string(REPLACE "<" "-" filename ${filename})
    string(REPLACE ">" "-" filename ${filename})
    string(REPLACE "," "-" filename ${filename})
    string(REPLACE " " ""  filename ${filename})
    message( " -- Generating file: ${filename}")
    configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/MultiFluidUpdate.cpp.template
                    ${filename} )
    list( APPEND physicsSolvers_sources ${filename} )                    

    foreach( NUM_COMP 1 2 3 4 5 )
        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/DirichletFaceBasedAssemblyKernel_${NUM_COMP}_${FLUID_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " ""  filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/DirichletFaceBasedAssemblyKernel.cpp.template
                        ${filename} )
        list( APPEND physicsSolvers_sources ${filename} )
    endforeach()
endforeach()

