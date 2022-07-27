/*
 * PVTDriverRunTestDeadOilFluid.cpp
 *
 *  Created on: Jul 27, 2022
 *      Author: settgast
 */



#include "PVTDriverRunTest.hpp"
#include "DeadOilFluid.hpp"
#include "BlackOilFluid.hpp"


namespace geosx
{
template void PVTDriver::runTest<constitutive::DeadOilFluid>( constitutive::DeadOilFluid &, arrayView2d< real64 > const & );
template void PVTDriver::runTest<constitutive::BlackOilFluid>( constitutive::BlackOilFluid &, arrayView2d< real64 > const & );
}
