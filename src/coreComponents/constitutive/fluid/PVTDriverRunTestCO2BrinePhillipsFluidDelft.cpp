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
template void PVTDriver::runTest< constitutive::CO2BrinePhillipsFluidDelft >( constitutive::CO2BrinePhillipsFluidDelft &, arrayView2d< real64 > const & );
}
