//
// Created by root on 10/25/22.
//

#include "RelpermDriverRunTest.hpp"
#include "BrooksCoreyBakerRelativePermeability.hpp"


namespace geosx
{
    template void RelpermDriver::runTest< geosx::constitutive::BrooksCoreyBakerRelativePermeability >( geosx::constitutive::BrooksCoreyBakerRelativePermeability &, arrayView3d< real64 > const & );
}
