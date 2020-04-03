/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseBaseKernels.cpp
 */

#include "SinglePhaseBaseKernels.hpp"

namespace geosx
{

namespace SinglePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

void
MobilityKernel::Compute( real64 const & dens,
                         real64 const & dDens_dPres,
                         real64 const & visc,
                         real64 const & dVisc_dPres,
                         real64 & mob,
                         real64 & dMob_dPres )
{
  mob = dens / visc;
  dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
}

void
MobilityKernel::Compute( real64 const & dens,
                         real64 const & visc,
                         real64 & mob )
{
  mob = dens / visc;
}

void MobilityKernel::Launch( localIndex const size,
                             arrayView2d< real64 const > const & dens,
                             arrayView2d< real64 const > const & dDens_dPres,
                             arrayView2d< real64 const > const & visc,
                             arrayView2d< real64 const > const & dVisc_dPres,
                             arrayView1d< real64 > const & mob,
                             arrayView1d< real64 > const & dMob_dPres )
{
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( SortedArrayView< localIndex const > targetSet,
                             arrayView2d< real64 const > const & dens,
                             arrayView2d< real64 const > const & dDens_dPres,
                             arrayView2d< real64 const > const & visc,
                             arrayView2d< real64 const > const & dVisc_dPres,
                             arrayView1d< real64 > const & mob,
                             arrayView1d< real64 > const & dMob_dPres )
{
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( localIndex const size,
                             arrayView2d< real64 const > const & dens,
                             arrayView2d< real64 const > const & visc,
                             arrayView1d< real64 > const & mob )
{
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

void MobilityKernel::Launch( SortedArrayView< localIndex const > targetSet,
                             arrayView2d< real64 const > const & dens,
                             arrayView2d< real64 const > const & visc,
                             arrayView1d< real64 > const & mob )
{
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

} // namespace SinglePhaseBaseKernels

} // namespace geosx
