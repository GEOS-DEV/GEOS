/*
 * PVTDriverRunTestExampleFluidModel.cpp
 *
 *  Created on: Jul 27, 2022
 *      Author: settgast
 */



#include "PVTDriverRunTest.hpp"
#include "ExampleFluidModel.hpp"


namespace geosx
{
template void PVTDriver::runTest< constitutive::ExampleFluidModel >( constitutive::ExampleFluidModel &, arrayView2d< real64 > const & );
}
