/*
 * PVTDriverRunTestDeadOilFluid.cpp
 *
 *  Created on: Jul 27, 2022
 *      Author: settgast
 */



#include "PVTDriverRunTest.hpp"
#include "CompositionalMultiphaseFluid.hpp"


namespace geosx
{
template void PVTDriver::runTest<constitutive::CompositionalMultiphaseFluid>( constitutive::CompositionalMultiphaseFluid &, arrayView2d< real64 > const & );
}
