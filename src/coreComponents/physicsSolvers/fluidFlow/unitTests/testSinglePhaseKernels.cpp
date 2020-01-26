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
#include "common/Logger.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseKernels.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testFlowKernelHelpers.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::SinglePhaseKernels;


TEST( SinglePhaseKernels, mobility )
{
  int constexpr NTEST = 3;

  real64 const dens[NTEST]        = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST] = { 1e-5,  1e-10,  0.0    };
  real64 const visc[NTEST]        = { 5.0,   2.0,    1.0    };
  real64 const dVisc_dPres[NTEST] = { 1e-7,  0.0,    0.0    };

  for (int i = 0; i < NTEST; ++i)
  {
    SCOPED_TRACE( "Input # " + std::to_string(i) );

    real64 mob;
    real64 dMob_dPres;

    MobilityKernel::Compute( dens[i], dDens_dPres[i], visc[i], dVisc_dPres[i], mob, dMob_dPres );

    // compute etalon
    real64 const mob_et = dens[i] / visc[i];
    real64 const dMob_dPres_et = mob_et * (dDens_dPres[i] / dens[i] - dVisc_dPres[i] / visc[i]);

    EXPECT_DOUBLE_EQ( mob, mob_et );
    EXPECT_DOUBLE_EQ( dMob_dPres, dMob_dPres_et );
  }
}

/**
 * @brief Test the assembly of accumulation contribution
 *
 * @note This only tests uncoupled version.
 * In future, porosity update will be elsewhere and this will be simplified.
 */
TEST( SinglePhaseKernels, accumulation )
{
  int constexpr NTEST = 3;

  real64 const densOld[NTEST]       = { 700.0, 1200.0, 1500.0 };
  real64 const densNew[NTEST]       = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST]   = { 1e-5,  1e-10,  0.0    };
  real64 const dVol[NTEST]          = { 1.0,   2.0,    10.0   };
  real64 const poroOld[NTEST]       = { 6e-2,  1e-1,   2e-1   };
  real64 const poroRef[NTEST]       = { 5e-2,  2e-1,   3e-1   };
  real64 const pvMult[NTEST]        = { 1.20,  1.10,   1.05   };
  real64 const dPvMult_dPres[NTEST] = { 1e-5,  1e-7,   0.0    };
  real64 const volume               = 1.0;


  for (int i = 0; i < NTEST; ++i)
  {
    SCOPED_TRACE( "Input # " + std::to_string(i) );

    real64 accum;
    real64 accumJacobian;
    real64 poroNew;

    AccumulationKernel<CellElementSubRegion>::Compute<false>( 0.0, densNew[i], densOld[i], dDens_dPres[i], volume, dVol[i],
                                        poroRef[i], poroOld[i], pvMult[i], dPvMult_dPres[i],
                                        0.0, 0.0, 0.0, 0.0, poroNew, accum, accumJacobian );

    // compute etalon
    real64 const poroNew_et = poroRef[i] * pvMult[i];
    real64 const accum_et = poroNew_et * densNew[i] * (volume + dVol[i]) - poroOld[i] * densOld[i] * volume;
    real64 const accumJacobian_et = dPvMult_dPres[i] * poroRef[i] * densNew[i] * (volume + dVol[i])
                                  + poroNew_et * dDens_dPres[i] * (volume + dVol[i]);

    EXPECT_DOUBLE_EQ( poroNew, poroNew_et );
    EXPECT_DOUBLE_EQ( accum, accum_et );
    EXPECT_DOUBLE_EQ( accumJacobian, accumJacobian_et );
  }
}

template<localIndex stencilSize>
void computeFlux( arraySlice1d<real64 const> const & weight,
                  real64 const * pres,
                  real64 const * dPres,
                  real64 const * gravCoef,
                  real64 const * mob,
                  real64 const * dMob_dPres,
                  real64 const * dens,
                  real64 const * dDens_dPres,
                  real64 const dt,
                  real64 & flux,
                  real64 (& dFlux_dP)[stencilSize] )
{
  localIndex constexpr numElems = 2;

  real64 densMean = 0.0;
  real64 dDensMean_dP[stencilSize] {};
  for (localIndex i = 0; i < numElems; ++i)
  {
    densMean += 0.5 * dens[i];
    dDensMean_dP[i] = 0.5 * dDens_dPres[i];
  }
  real64 potDif = 0.0;
  real64 sumWeightGrav = 0;
  for (localIndex i = 0; i < stencilSize; ++i)
  {
    potDif += weight[i] * (pres[i] + dPres[i] - densMean * gravCoef[i]);
    sumWeightGrav += weight[i] * gravCoef[i];
  }
  localIndex const k_up = (potDif >= 0) ? 0 : 1;
  flux = dt * potDif * mob[k_up];
  for (localIndex i = 0; i < stencilSize; ++i)
  {
    dFlux_dP[i] = dt * ( weight[i] - sumWeightGrav * dDensMean_dP[i] ) * mob[k_up];
  }
  dFlux_dP[k_up] += dt * potDif * dMob_dPres[k_up];
}

template<bool FULL, localIndex stencilSize>
void testFluxKernel( CellElementStencilTPFA const & stencil,
                     real64 const * pres,
                     real64 const * dPres,
                     real64 const * gravCoef,
                     real64 const * mob,
                     real64 const * dMob_dPres,
                     real64 const * dens,
                     real64 const * dDens_dPres,
                     real64 const dt )
{
  localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  localIndex constexpr fluidIndex = 0;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename CellElementStencilTPFA::WeightContainerViewConstType const & weights = stencil.getWeights();

  auto presView        = AccessorHelper<FULL>::template makeElementAccessor<1> ( pres,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0] );
  auto dPresView       = AccessorHelper<FULL>::template makeElementAccessor<1> ( dPres,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0]);
  auto gravCoefView    = AccessorHelper<FULL>::template makeElementAccessor<1> ( gravCoef,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0]);
  auto mobView         = AccessorHelper<FULL>::template makeElementAccessor<1> ( mob,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0] );
  auto dMob_dPresView  = AccessorHelper<FULL>::template makeElementAccessor<1> ( dMob_dPres,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0] );
  auto densView        = AccessorHelper<FULL>::template makeMaterialAccessor<2>( dens,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0],
                                                                                 fluidIndex );
  auto dDens_dPresView = AccessorHelper<FULL>::template makeMaterialAccessor<2>( dDens_dPres,
                                                                                 stencilSize,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0],
                                                                                 fluidIndex );

  array1d<real64> flux( numElems );
  array2d<real64> fluxJacobian( numElems, stencilSize );



  FluxKernel::Compute( stencilSize,
                       seri[0],
                       sesri[0],
                       sei[0],
                       weights[0],
                       presView.toViewConst(),
                       dPresView.toViewConst(),
                       gravCoefView.toViewConst(),
                       densView.toViewConst(),
                       dDens_dPresView.toViewConst(),
                       mobView.toViewConst(),
                       dMob_dPresView.toViewConst(),
                       fluidIndex,
                       dt,
                       flux,
                       fluxJacobian );

  real64 flux_et;
  real64 dFlux_dP_et[stencilSize];

  // compute etalon
  computeFlux( weights[0],
               pres,
               dPres,
               gravCoef,
               mob,
               dMob_dPres,
               dens,
               dDens_dPres,
               dt,
               flux_et,
               dFlux_dP_et );

  EXPECT_DOUBLE_EQ( flux[0],  flux_et );
  EXPECT_DOUBLE_EQ( flux[1], -flux_et );
  for (localIndex i = 0; i < stencilSize; ++i)
  {
    EXPECT_DOUBLE_EQ( fluxJacobian[0][i],  dFlux_dP_et[i] );
    EXPECT_DOUBLE_EQ( fluxJacobian[1][i], -dFlux_dP_et[i] );
  }
}

TEST( SinglePhaseKernels, fluxFull )
{
  localIndex constexpr stencilSize = 2;


  CellElementStencilTPFA stencil;

  localIndex elemReg[2] = {0,1};
  localIndex elemSubReg[2] = {0,0};
  localIndex elemIndex[2] = {1,0};
  real64 weight[] = { 1e-12, -1e-12 };
  stencil.add( stencilSize,
               elemReg,
               elemSubReg,
               elemIndex,
               weight,
               0 );


  int constexpr NTEST = 3;

  // we keep these around for easy aggregate initialization
  real64  const presData       [NTEST][stencilSize] = { { 1e+6, 2e+6 }, { 2e+6, 2e+6 }, { 2e+6, 2e+6 } };
  real64  const dPresData      [NTEST][stencilSize] = { { 1e+5, 1e+5 }, { 1e+5, 2e+5 }, { 1e+5, 1e+5 } };
  real64  const gravCoefData   [NTEST][stencilSize] = { { 1e+3, 5e+2 }, { 1e+3, 1e+3 }, { 0.0,  1e+3 } };
  real64  const mobData        [NTEST][stencilSize] = { { 1e+6, 2e+6 }, { 2e+6, 1e+6 }, { 2e+6, 5e+6 } };
  real64  const dMob_dPresData [NTEST][stencilSize] = { { 1e-6, 2e-6 }, { 1e-6, 2e-6 }, { 1e-6, 2e-6 } };
  real64  const densData       [NTEST][stencilSize] = { { 1e+3, 2e+3 }, { 2e+3, 3e+3 }, { 2e+3, 1e+3 } };
  real64  const dDens_dPresData[NTEST][stencilSize] = { { 1e-6, 2e-6 }, { 2e-6, 3e-6 }, { 2e-6, 2e-6 } };
  real64  const dt             [NTEST]              = { 1.0,            1e+5,           1e+8           };


  for (int i = 0; i < NTEST; ++i)
  {
    SCOPED_TRACE( "Input # " + std::to_string(i) );

    testFluxKernel<true,2>( stencil,
                          presData[i],
                          dPresData[i],
                          gravCoefData[i],
                          mobData[i],
                          dMob_dPresData[i],
                          densData[i],
                          dDens_dPresData[i],
                          dt[i] );

  }
}

TEST( SinglePhaseKernels, fluxRegion )
{
  localIndex constexpr stencilSize = 2;
  CellElementStencilTPFA stencil;

  localIndex elemReg[2] = {0,0};
  localIndex elemSubReg[2] = {0,0};
  localIndex elemIndex[2] = {1,0};
  real64 weight[] = { 1e-12, -1e-12 };
  stencil.add( stencilSize,
               elemReg,
               elemSubReg,
               elemIndex,
               weight,
               0 );

  int constexpr NTEST = 3;

  // we keep these around for easy aggregate initialization
  real64  const presData       [NTEST][stencilSize] = { { 1e+6, 2e+6 }, { 2e+6, 2e+6 }, { 2e+6, 2e+6 } };
  real64  const dPresData      [NTEST][stencilSize] = { { 1e+5, 1e+5 }, { 1e+5, 2e+5 }, { 1e+5, 1e+5 } };
  real64  const gravCoefData   [NTEST][stencilSize] = { { 1e+3, 5e+2 }, { 1e+3, 1e+3 }, { 0.0,  1e+3 } };
  real64  const mobData        [NTEST][stencilSize] = { { 1e+6, 2e+6 }, { 2e+6, 1e+6 }, { 2e+6, 5e+6 } };
  real64  const dMob_dPresData [NTEST][stencilSize] = { { 1e-6, 2e-6 }, { 1e-6, 2e-6 }, { 1e-6, 2e-6 } };
  real64  const densData       [NTEST][stencilSize] = { { 1e+3, 2e+3 }, { 2e+3, 3e+3 }, { 2e+3, 1e+3 } };
  real64  const dDens_dPresData[NTEST][stencilSize] = { { 1e-6, 2e-6 }, { 2e-6, 3e-6 }, { 2e-6, 2e-6 } };
  real64  const dt             [NTEST]              = { 1.0,            1e+5,           1e+8           };


  for (int i = 0; i < NTEST; ++i)
  {
    SCOPED_TRACE( "Input # " + std::to_string(i) );

    testFluxKernel<false,2>( stencil,
                           presData[i],
                           dPresData[i],
                           gravCoefData[i],
                           mobData[i],
                           dMob_dPresData[i],
                           densData[i],
                           dDens_dPresData[i],
                           dt[i] );

  }
}

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}

