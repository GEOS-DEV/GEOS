/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MobilityKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_MOBILITYKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_MOBILITYKERNEL_HPP


namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  GEOS_HOST_DEVICE
  inline
  static void
  compute( geos::real64 const & dens,
           geos::real64 const & dDens_dPres,
           geos::real64 const & dDens_dTemp,
           geos::real64 const & visc,
           geos::real64 const & dVisc_dPres,
           geos::real64 const & dVisc_dTemp,
           geos::real64 & mob,
           geos::real64 & dMob_dPres,
           geos::real64 & dMob_dTemp )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
    dMob_dTemp = dDens_dTemp / visc - mob / visc * dVisc_dTemp;
  }

  GEOS_HOST_DEVICE
  inline
  static void
  compute( geos::real64 const & dens,
           geos::real64 const & visc,
           geos::real64 & mob )
  {
    mob = dens / visc;
  }

  template< typename POLICY >
  static void launch( geos::localIndex const size,
                      geos::arrayView2d< geos::real64 const > const & dens,
                      geos::arrayView2d< geos::real64 const > const & dDens_dPres,
                      geos::arrayView2d< geos::real64 const > const & dDens_dTemp,
                      geos::arrayView2d< geos::real64 const > const & visc,
                      geos::arrayView2d< geos::real64 const > const & dVisc_dPres,
                      geos::arrayView2d< geos::real64 const > const & dVisc_dTemp,
                      geos::arrayView1d< geos::real64 > const & mob,
                      geos::arrayView1d< geos::real64 > const & dMob_dPres,
                      geos::arrayView1d< geos::real64 > const & dMob_dTemp )
  {
    geos::forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( geos::localIndex const a )
    {
      compute( dens[a][0],
               dDens_dPres[a][0],
               dDens_dTemp[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               dVisc_dTemp[a][0],
               mob[a],
               dMob_dPres[a],
               dMob_dTemp[a] );
    } );
  }

  template< typename POLICY >
  static void launch( geos::localIndex const size,
                      geos::arrayView2d< geos::real64 const > const & dens,
                      geos::arrayView2d< geos::real64 const > const & visc,
                      geos::arrayView1d< geos::real64 > const & mob )
  {
    geos::forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( geos::localIndex const a )
    {
      compute( dens[a][0],
               visc[a][0],
               mob[a] );
    } );
  }
};

}

}

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_MOBILITYKERNEL_HPP
