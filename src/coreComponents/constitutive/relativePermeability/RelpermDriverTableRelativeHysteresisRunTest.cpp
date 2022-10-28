//
// Created by root on 10/25/22.
//

#include "RelpermDriverRunTest.hpp"
#include "TableRelativePermeabilityHysteresis.hpp"


namespace geosx
{
    template void RelpermDriver::runTest< geosx::constitutive::TableRelativePermeabilityHysteresis >( geosx::constitutive::TableRelativePermeabilityHysteresis &, arrayView3d< real64 > const & );
}
