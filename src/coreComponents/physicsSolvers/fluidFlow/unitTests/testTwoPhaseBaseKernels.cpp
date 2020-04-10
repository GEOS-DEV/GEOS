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

// Source includes
#include "managers/initialization.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseBaseKernels.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::TwoPhaseBaseKernels;
using namespace geosx::testing;

static real64 constexpr relTol = 1e-8;

TEST( TwoPhaseBaseKernels, phaseMobility )
{
  int constexpr NTEST = 3;
  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  array2d< real64 > phaseDens, dPhaseDens_dPres;
  phaseDens.resize( NTEST, numPhases ); dPhaseDens_dPres.resize( NTEST, numPhases );
  phaseDens( 0, 0 ) = 800.0; phaseDens( 1, 0 ) = 600.0; phaseDens( 2, 0 ) = 10.0;
  phaseDens( 0, 1 ) = 8.0; phaseDens( 1, 1 ) = 1200.0; phaseDens( 2, 1 ) = 50.0;
  dPhaseDens_dPres( 0, 0 ) = 1e-5; dPhaseDens_dPres( 1, 0 ) = 1e-7; dPhaseDens_dPres( 2, 0 ) = 0.0;
  dPhaseDens_dPres( 0, 1 ) = 1e-8; dPhaseDens_dPres( 1, 1 ) = 0.0; dPhaseDens_dPres( 2, 1 ) = 1e-12;

  array2d< real64 > phaseVisc, dPhaseVisc_dPres;
  phaseVisc.resize( NTEST, numPhases ); dPhaseVisc_dPres.resize( NTEST, numPhases );
  phaseVisc( 0, 0 ) = 0.01; phaseVisc( 1, 0 ) = 0.001; phaseVisc( 2, 0 ) = 0.007;
  phaseVisc( 0, 1 ) = 0.006; phaseVisc( 1, 1 ) = 0.02; phaseVisc( 2, 1 ) = 0.0001;
  dPhaseVisc_dPres( 0, 0 ) = 1e-5; dPhaseVisc_dPres( 1, 0 ) = 1e-7; dPhaseVisc_dPres( 2, 0 ) = 0.0;
  dPhaseVisc_dPres( 0, 1 ) = 1e-8; dPhaseVisc_dPres( 1, 1 ) = 0.0; dPhaseVisc_dPres( 2, 1 ) = 1e-12;

  array2d< real64 > phaseRelPerm;
  array3d< real64 > dPhaseRelPerm_dSat;
  phaseRelPerm.resize( NTEST, numPhases );
  dPhaseRelPerm_dSat.resizeDimension< 0 >( NTEST );
  dPhaseRelPerm_dSat.resizeDimension< 1 >( numPhases );
  dPhaseRelPerm_dSat.resizeDimension< 2 >( numPhases );
  phaseRelPerm( 0, 0 ) = 0.45; phaseRelPerm( 1, 0 ) = 0.3; phaseRelPerm( 2, 0 ) = 0.7;
  phaseRelPerm( 0, 1 ) = 0.6; phaseRelPerm( 1, 1 ) = 0.9; phaseRelPerm( 2, 1 ) = 0.4;
  dPhaseRelPerm_dSat[0][0][0] = 0.25; dPhaseRelPerm_dSat[1][0][0] = 1e-7; dPhaseRelPerm_dSat[2][0][0] = 0.8;
  dPhaseRelPerm_dSat[0][0][1] = 0.0; dPhaseRelPerm_dSat[1][0][1] = 0; dPhaseRelPerm_dSat[2][0][1] = 0;
  dPhaseRelPerm_dSat[0][1][0] =  0.0; dPhaseRelPerm_dSat[1][1][0] = 0.0; dPhaseRelPerm_dSat[2][1][0] = 0.0;
  dPhaseRelPerm_dSat[0][1][1] = -0.5; dPhaseRelPerm_dSat[1][1][1] = -0.75; dPhaseRelPerm_dSat[2][1][1] = -0.2;

  array1d< real64 > phaseMobility, dPhaseMobility_dPres, dPhaseMobility_dSat;
  phaseMobility.resize( numPhases );
  dPhaseMobility_dPres.resize( numPhases );
  dPhaseMobility_dSat.resize( numPhases );

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    PhaseMobilityKernel::Compute( phaseDens[i],
                                  dPhaseDens_dPres[i],
                                  phaseVisc[i],
                                  dPhaseVisc_dPres[i],
                                  phaseRelPerm[i],
                                  dPhaseRelPerm_dSat[i],
                                  phaseMobility,
                                  dPhaseMobility_dPres,
                                  dPhaseMobility_dSat );

    // compute etalon
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const relPermOverVisc = phaseRelPerm( i, ip ) / phaseVisc( i, ip );
      real64 const dRelPermOverVisc_dPres = -dPhaseVisc_dPres( i, ip ) * phaseRelPerm( i, ip )
                                            / ( phaseVisc( i, ip ) * phaseVisc( i, ip ) );
      real64 const dRelPermOverVisc_dSat = dPhaseRelPerm_dSat[i][ip][ip] / phaseVisc( i, ip );

      real64 const phaseMobility_et = phaseDens( i, ip ) * relPermOverVisc;
      real64 const dPhaseMobility_dPres_et = dPhaseDens_dPres( i, ip ) * relPermOverVisc
                                             + phaseDens( i, ip ) * dRelPermOverVisc_dPres;
      real64 dPhaseMobility_dSat_et = phaseDens( i, ip ) * dRelPermOverVisc_dSat;

      if( ip == 1 )
      {
        dPhaseMobility_dSat_et *= -1;
      }

      EXPECT_DOUBLE_EQ( phaseMobility( ip ), phaseMobility_et );
      EXPECT_DOUBLE_EQ( dPhaseMobility_dPres( ip ), dPhaseMobility_dPres_et );
      EXPECT_DOUBLE_EQ( dPhaseMobility_dSat( ip ), dPhaseMobility_dSat_et );
    }
  }
}

TEST( TwoPhaseBaseKernels, accumulation )
{
  int constexpr NTEST = 3;

  localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

  array2d< real64 > phaseDensNew, phaseDensOld, dPhaseDensNew_dPres;
  phaseDensNew.resize( NTEST, numPhases ); phaseDensOld.resize( NTEST, numPhases ); dPhaseDensNew_dPres.resize( NTEST, numPhases );
  phaseDensNew( 0, 0 ) = 800.0; phaseDensNew( 1, 0 ) = 600.0; phaseDensNew( 2, 0 ) = 10.0;
  phaseDensNew( 0, 1 ) = 8.0; phaseDensNew( 1, 1 ) = 1200.0; phaseDensNew( 2, 1 ) = 50.0;
  phaseDensOld( 0, 0 ) = 843.0; phaseDensOld( 1, 0 ) = 550.0; phaseDensOld( 2, 0 ) = 20.0;
  phaseDensOld( 0, 1 ) = 11.0; phaseDensOld( 1, 1 ) = 1100.0; phaseDensOld( 2, 1 ) = 100.0;
  dPhaseDensNew_dPres( 0, 0 ) = 1e-5; dPhaseDensNew_dPres( 1, 0 ) = 1e-7; dPhaseDensNew_dPres( 2, 0 ) = 0.0;
  dPhaseDensNew_dPres( 0, 1 ) = 1e-8; dPhaseDensNew_dPres( 1, 1 ) = 0.0; dPhaseDensNew_dPres( 2, 1 ) = 1e-12;

  array2d< real64 > phaseSat, dPhaseSat;
  phaseSat.resize( NTEST, numPhases ); dPhaseSat.resize( NTEST, numPhases );
  phaseSat( 0, 0 ) = 0.1; phaseSat( 1, 0 ) = 0.9; phaseSat( 2, 0 ) = 0.5;
  phaseSat( 0, 1 ) = 0.8; phaseSat( 1, 1 ) = 0.05; phaseSat( 2, 1 ) = 0.3;
  dPhaseSat( 0, 0 ) = 0.05; dPhaseSat( 1, 0 ) = 0.01; dPhaseSat( 2, 0 ) = 0.05;
  dPhaseSat( 0, 1 ) = 0.05; dPhaseSat( 1, 1 ) = 0.04; dPhaseSat( 2, 1 ) = 0.15;

  real64 const dVol[NTEST]          = { 0.0, 0.0, 0.0   };
  real64 const poroOld[NTEST]       = { 6e-2, 1e-1, 2e-1   };
  real64 const poroRef[NTEST]       = { 5e-2, 2e-1, 3e-1   };
  real64 const pvMult[NTEST]        = { 1.20, 1.10, 1.05   };
  real64 const dPvMult_dPres[NTEST] = { 1e-5, 1e-7, 0.0    };
  real64 const volume               = 1.0;

  array1d< real64 > accum;
  array2d< real64 > accumJacobian;
  accum.resize( numPhases );
  accumJacobian.resize( numPhases, numPhases );

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    AccumulationKernel::Compute( volume + dVol[i],
                                 poroOld[i],
                                 poroRef[i],
                                 pvMult[i],
                                 dPvMult_dPres[i],
                                 phaseSat[i],
                                 dPhaseSat[i],
                                 phaseDensOld[i],
                                 phaseDensNew[i],
                                 dPhaseDensNew_dPres[i],
                                 accum,
                                 accumJacobian );

    // compute etalon
    real64 const poroNew_et = poroRef[i] * pvMult[i];

    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const phaseSatNew_et = phaseSat( i, ip ) + dPhaseSat( i, ip );
      real64 const dPhaseSatNew_et = (ip == 0) ? 1 : -1;
      real64 const accum_et = poroNew_et * phaseSatNew_et * phaseDensNew( i, ip ) * (volume + dVol[i])
                              - poroOld[i] * phaseSat( i, ip ) * phaseDensOld( i, ip ) * volume;
      real64 const dAccum_dPres_et = dPvMult_dPres[i] * poroRef[i] * phaseSatNew_et * phaseDensNew( i, ip ) * (volume + dVol[i])
                                     + poroNew_et * phaseSatNew_et * dPhaseDensNew_dPres( i, ip ) * (volume + dVol[i]);
      real64 const dAccum_dSat_et = poroNew_et * dPhaseSatNew_et * phaseDensNew( i, ip ) * (volume + dVol[i]);

      checkRelativeError( accum[ip], accum_et, relTol );
      checkRelativeError( accumJacobian[ip][0], dAccum_dPres_et, relTol );
      checkRelativeError( accumJacobian[ip][1], dAccum_dSat_et, relTol );
    }
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
