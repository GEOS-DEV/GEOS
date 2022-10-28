//
// Created by root on 10/25/22.
//

#include "RelpermDriverRunTest.hpp"
#include "VanGenuchtenBakerRelativePermeability.hpp"


namespace geosx
{
    template void RelpermDriver::runTest< geosx::constitutive::VanGenuchtenBakerRelativePermeability >( geosx::constitutive::VanGenuchtenBakerRelativePermeability &, arrayView3d< real64 > const & );
}
