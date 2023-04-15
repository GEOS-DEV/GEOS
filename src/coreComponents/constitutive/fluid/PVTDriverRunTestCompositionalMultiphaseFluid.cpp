/*
 * PVTDriverRunTestDeadOilFluid.cpp
 *
 *  Created on: Jul 27, 2022
 *      Author: settgast
 */



#include "PVTDriverRunTest.hpp"
#include "CompositionalMultiphaseFluid.hpp"


namespace geos
{
template void PVTDriver::runTest< constitutive::CompositionalMultiphaseFluid >( constitutive::CompositionalMultiphaseFluid &, arrayView2d< real64 > const & );
}
