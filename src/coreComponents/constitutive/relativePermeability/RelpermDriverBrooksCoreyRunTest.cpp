//
// Created by root on 10/25/22.
//

#include "RelpermDriverRunTest.hpp"
#include "BrooksCoreyRelativePermeability.hpp"


namespace geosx
{
    template void RelpermDriver::runTest< geosx::constitutive::BrooksCoreyRelativePermeability >( geosx::constitutive::BrooksCoreyRelativePermeability &, arrayView3d< real64 > const & );
}
