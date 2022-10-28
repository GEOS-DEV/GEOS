//
// Created by root on 10/25/22.
//

#include "RelpermDriverRunTest.hpp"
#include "TableRelativePermeability.hpp"


namespace geosx
{
    template void RelpermDriver::runTest< geosx::constitutive::TableRelativePermeability >( geosx::constitutive::TableRelativePermeability &, arrayView3d< real64 > const & );
}
