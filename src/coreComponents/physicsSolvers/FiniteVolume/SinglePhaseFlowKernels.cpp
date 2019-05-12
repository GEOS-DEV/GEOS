/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SinglePhaseFlowKernels.cpp
 */

#include "SinglePhaseFlowKernels.hpp"

namespace geosx
{

namespace SinglePhaseFlowKernels
{

/******************************** MobilityKernel ********************************/

inline RAJA_HOST_DEVICE void
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

inline RAJA_HOST_DEVICE void
MobilityKernel::Compute( real64 const & dens,
                         real64 const & visc,
                         real64 & mob )
{
  mob = dens / visc;
}

void MobilityKernel::Launch( localIndex begin, localIndex end,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & dDens_dPres,
                             arrayView2d<real64 const> const & visc,
                             arrayView2d<real64 const> const & dVisc_dPres,
                             arrayView1d<real64> const & mob,
                             arrayView1d<real64> const & dMob_dPres )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( set<localIndex> targetSet,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & dDens_dPres,
                             arrayView2d<real64 const> const & visc,
                             arrayView2d<real64 const> const & dVisc_dPres,
                             arrayView1d<real64> const & mob,
                             arrayView1d<real64> const & dMob_dPres )
{
  forall_in_set( targetSet.values(), targetSet.size(), GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( localIndex begin, localIndex end,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & visc,
                             arrayView1d<real64> const & mob )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

void MobilityKernel::Launch( set<localIndex> targetSet,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & visc,
                             arrayView1d<real64> const & mob )
{
  forall_in_set( targetSet.values(), targetSet.size(), GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

} // namespace SinglePhaseFlowKernels

} // namespace geosx
