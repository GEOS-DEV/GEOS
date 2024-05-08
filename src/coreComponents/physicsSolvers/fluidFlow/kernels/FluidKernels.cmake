set( kernelPath "coreComponents/physicsSolvers/fluidFlow/kernels" )
set( fluidTypes DeadOilFluid
                BlackOilFluid
                CompositionalTwoPhasePengRobinsonConstantViscosity
                CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity
                CompositionalTwoPhasePengRobinsonLBCViscosity
                CompositionalTwoPhaseSoaveRedlichKwongLBCViscosity
                CO2BrinePhillipsFluid
                CO2BrineEzrokhiFluid
                CO2BrinePhillipsThermalFluid
                CO2BrineEzrokhiThermalFluid  )

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
endforeach()

