/*
 * PVTDriverRunTestDeadOilFluid.cpp
 *
 *  Created on: Jul 27, 2022
 *      Author: settgast
 */



#include "PVTDriverRunTest.hpp"
#include "CO2BrineFluid.hpp"


namespace geosx
{
template void PVTDriver::runTest<constitutive::CO2BrinePhillipsFluid>( constitutive::CO2BrinePhillipsFluid &, arrayView2d< real64 > const & );
template void PVTDriver::runTest<constitutive::CO2BrineEzrokhiFluid>( constitutive::CO2BrineEzrokhiFluid &, arrayView2d< real64 > const & );
template void PVTDriver::runTest<constitutive::CO2BrinePhillipsThermalFluid>( constitutive::CO2BrinePhillipsThermalFluid &, arrayView2d< real64 > const & );
template void PVTDriver::runTest<constitutive::CO2BrineEzrokhiThermalFluid>( constitutive::CO2BrineEzrokhiThermalFluid &, arrayView2d< real64 > const & );
}
