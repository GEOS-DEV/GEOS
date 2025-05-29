/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

/**
 * This test compares the solubility calculated from our implementation of the Soreide Whitson
 * model against published experimental data.
 */
using SolubilityData = std::tuple<
  real64 const,     // pressure
  real64 const,     // temperature
  real64 const,     // salinity
  integer const,    // component (see TestFluid.hpp)
  real64 const,     // measured (experimental) liquid mole fraction
  real64 const      // model liquid mole fraction (our implementation)
  >;

class SoreideWhitsonSolubilityTestFixture : public ::testing::TestWithParam< SolubilityData >
{
protected:
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr real64 expRelTol = 0.2;  // We can only match all the experimental data to 20% accurace
  static constexpr real64 expAbsTol = 1.0;  // Ignore absolute differences
};

TEST_P( SoreideWhitsonSolubilityTestFixture, testSolubility )
{
  auto const [pressure, temperature, salinity, component, measuredMoleFraction, expectedMoleFraction] = GetParam();

  constitutive::compositional::FlashData flashData;
  flashData.liquidEos = EquationOfStateType::SoreideWhitson;
  flashData.vapourEos = EquationOfStateType::PengRobinson;
  flashData.salinity = salinity;

  // Create a 2-component fluid
  static constexpr integer numComps = 2;
  auto fluid = TestFluid< numComps >::create( {component, Fluid::H2O} );

  /* UNCRUSTIFY-OFF */
  // Soreide-Whitson correlations work only with "true" values of component parameters
  // kij_NA is the binary interation coefficient in the gas phase (see Table 5 Soreide-Whitson (1992))
  std::unordered_map<integer, std::array<real64 const, 6> const> const componentDatabase = {
    //             Mw            Pc            Tc            Vc            Ac           kij_NA
    {Fluid::H2O, { 1.80153e-02,  2.20640e+07,  6.47096e+02,  5.59480e-05,  3.44300e-01, 0.0000 }},
    {Fluid::CO2, { 4.40095e-02,  7.37730e+06,  3.04128e+02,  9.41185e-05,  2.23940e-01, 0.1896 }},
    {Fluid::N2,  { 2.80134e-02,  3.39580e+06,  1.26192e+02,  8.94142e-05,  3.72000e-02, 0.4778 }},
    {Fluid::H2S, { 3.40809e-02,  9.00000e+06,  3.73100e+02,  9.81354e-05,  1.00500e-01, 0.0000 }},
    {Fluid::C1,  { 1.60425e-02,  4.59920e+06,  1.90564e+02,  9.86278e-05,  1.14200e-02, 0.4850 }},
    {Fluid::C2,  { 3.00690e-02,  4.87220e+06,  3.05322e+02,  1.45839e-04,  9.95000e-02, 0.4920 }},
    {Fluid::H2,  { 2.01588e-03,  1.29640e+06,  3.31450e+01,  6.44828e-05, -2.19000e-01, 0.0000 }},
  };
  /* UNCRUSTIFY-ON */

  auto const & componentData = componentDatabase.find( component )->second;
  fluid->molecularWeight[0] = componentData[0];
  fluid->criticalPressure[0] = componentData[1];
  fluid->criticalTemperature[0] = componentData[2];
  fluid->criticalVolume[0] = componentData[3];
  fluid->acentricFactor[0] = componentData[4];
  auto const & h2oData = componentDatabase.find( Fluid::H2O )->second;
  fluid->molecularWeight[1] = h2oData[0];
  fluid->criticalPressure[1] = h2oData[1];
  fluid->criticalTemperature[1] = h2oData[2];
  fluid->criticalVolume[1] = h2oData[3];
  fluid->acentricFactor[1] = h2oData[4];

  real64 kij_NA = componentData[5];     // Vapour binary interaction coefficient
  real64 max_x_c = 0.08;
  if( component == Fluid::CO2 )
  {
    max_x_c = 0.2;     // CO2 is more soluble
  }
  else if( component == Fluid::H2S )
  {
    real64 const Tr = temperature / fluid->criticalTemperature[0];
    kij_NA = 0.19031 - 0.05965*Tr;     /* Eqn 17. Soreide-Whitson (1992) */
  }

  fluid->setBinaryCoefficients( Feed< 1 >{kij_NA} );
  auto componentProperties = fluid->createKernelWrapper();

  real64 vapourFraction = -1.0;
  stackArray2d< real64, 3*numComps > compositions( 3, numComps );
  auto totalComposition = compositions[0];
  auto liquidComposition = compositions[1];
  auto vapourComposition = compositions[2];
  stackArray2d< real64, numComps > kValues( 1, numComps );

  totalComposition[0] = max_x_c;
  totalComposition[1] = 1.0 - totalComposition[0];

  kValues( 0, 0 ) = 100.0;
  kValues( 0, 1 ) = 0.01;

  bool status = NegativeTwoPhaseFlash::compute(
    numComps,
    pressure,
    temperature,
    totalComposition.toSliceConst(),
    componentProperties,
    flashData,
    kValues.toSlice(),
    vapourFraction,
    liquidComposition,
    vapourComposition );

  // Check the flash success result
  ASSERT_TRUE( status );

  if( !status )
  {
    return;
  }

  real64 const moleFraction = liquidComposition[0];
  checkRelativeError( measuredMoleFraction, moleFraction, expRelTol, expAbsTol );
  checkRelativeError( expectedMoleFraction, moleFraction, relTol, absTol );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */
// References:
// - Chabab et al. (2024) https://doi.org/10.1016/j.ijhydene.2023.10.290
// - Gillespie et al. (1984) Wilson, G. M., Gillespie, P. C., & Owens, J. L. (1985). Sour water
//       equilibria extended to high temperatures and with inerts present. In Annual convention
//       of Gas Processors Association. 64 (pp. 282-288).
// - Lee & Mather (1977) https://inis.iaea.org/records/wwxbh-09x36
// - O'Sullivan & Smith (1970) https://doi.org/10.1021/j100702a012
// - Selleck et al. (1952) https://doi.org/10.1021/ie50513a064
// - Soreide & Whitson (1992) https://doi.org/10.1016/0378-3812(92)85105-H
// - Takenouchi & Kennedy (1964) https://doi.org/10.2475/ajs.262.9.1055
INSTANTIATE_TEST_SUITE_P( SoreideWhitsonFlashTest, SoreideWhitsonSolubilityTestFixture,
  ::testing::ValuesIn< SolubilityData >({
    // CH4 data from O'Sullivan & Smith (1970)
    // Digitized from Fig 14 of Soreide & Whitson (1992)
    { 1.0039e+07, 376.15, 0.0, Fluid::C1,  1.352e-03, 1.246435e-03 }, 
    { 2.0079e+07, 376.15, 0.0, Fluid::C1,  2.197e-03, 2.139535e-03 }, 
    { 3.0118e+07, 376.15, 0.0, Fluid::C1,  2.873e-03, 2.815925e-03 }, 
    { 4.0157e+07, 376.15, 0.0, Fluid::C1,  3.320e-03, 3.359950e-03 }, 
    { 5.0059e+07, 376.15, 0.0, Fluid::C1,  3.847e-03, 3.810142e-03 }, 
    { 6.0373e+07, 376.15, 0.0, Fluid::C1,  4.195e-03, 4.213468e-03 }, 
    { 2.0079e+07, 376.15, 1.0, Fluid::C1,  1.690e-03, 1.667933e-03 }, 
    { 3.0118e+07, 376.15, 1.0, Fluid::C1,  2.207e-03, 2.190537e-03 }, 
    { 4.0295e+07, 376.15, 1.0, Fluid::C1,  2.565e-03, 2.614429e-03 }, 
    { 5.0196e+07, 376.15, 1.0, Fluid::C1,  2.893e-03, 2.958925e-03 }, 
    { 6.0373e+07, 376.15, 1.0, Fluid::C1,  3.201e-03, 3.262970e-03 }, 
    { 2.0216e+07, 376.15, 4.0, Fluid::C1,  8.350e-04, 7.711673e-04 }, 
    { 3.0118e+07, 376.15, 4.0, Fluid::C1,  1.083e-03, 1.003461e-03 }, 
    { 4.0157e+07, 376.15, 4.0, Fluid::C1,  1.223e-03, 1.190968e-03 }, 
    { 5.0196e+07, 376.15, 4.0, Fluid::C1,  1.332e-03, 1.346522e-03 }, 
    { 6.0236e+07, 376.15, 4.0, Fluid::C1,  1.451e-03, 1.479266e-03 }, 
    // CO2 data from Takenouchi & Kennedy (1964)
    // Digitized from Fig 16 of Soreide & Whitson (1992)
    { 1.0031e+07, 423.15, 0.0, Fluid::CO2, 1.290e-02, 1.177661e-02 }, 
    { 2.0179e+07, 423.15, 0.0, Fluid::CO2, 2.148e-02, 2.041963e-02 }, 
    { 3.0034e+07, 423.15, 0.0, Fluid::CO2, 2.580e-02, 2.623686e-02 }, 
    { 3.9947e+07, 423.15, 0.0, Fluid::CO2, 2.943e-02, 3.069806e-02 }, 
    { 1.0089e+07, 423.15, 1.1, Fluid::CO2, 1.176e-02, 1.012953e-02 }, 
    { 2.0061e+07, 423.15, 1.1, Fluid::CO2, 1.892e-02, 1.720937e-02 }, 
    { 3.0034e+07, 423.15, 1.1, Fluid::CO2, 2.239e-02, 2.202448e-02 }, 
    { 4.0006e+07, 423.15, 1.1, Fluid::CO2, 2.477e-02, 2.562682e-02 }, 
    { 1.0089e+07, 423.15, 4.3, Fluid::CO2, 4.205e-03, 4.186474e-03 }, 
    { 2.0003e+07, 423.15, 4.3, Fluid::CO2, 7.557e-03, 6.810234e-03 }, 
    { 2.9975e+07, 423.15, 4.3, Fluid::CO2, 1.040e-02, 8.483147e-03 }, 
    { 4.0064e+07, 423.15, 4.3, Fluid::CO2, 1.182e-02, 9.681339e-03 }, 
    // N2 data from O'Sullivan & Smith (1970)
    // Digitized from Fig 15 of Soreide & Whitson (1992)
    { 1.0139e+07, 376.15, 0.0, Fluid::N2,  7.621e-04, 8.008581e-04 }, 
    { 2.0279e+07, 376.15, 0.0, Fluid::N2,  1.425e-03, 1.487562e-03 }, 
    { 3.0320e+07, 376.15, 0.0, Fluid::N2,  1.994e-03, 2.079075e-03 }, 
    { 4.0265e+07, 376.15, 0.0, Fluid::N2,  2.507e-03, 2.598776e-03 }, 
    { 5.0306e+07, 376.15, 0.0, Fluid::N2,  2.963e-03, 3.070291e-03 }, 
    { 6.0348e+07, 376.15, 0.0, Fluid::N2,  3.362e-03, 3.497819e-03 }, 
    { 1.0139e+07, 376.15, 1.0, Fluid::N2,  5.912e-04, 5.746448e-04 }, 
    { 2.0181e+07, 376.15, 1.0, Fluid::N2,  1.090e-03, 1.060051e-03 }, 
    { 3.0320e+07, 376.15, 1.0, Fluid::N2,  1.524e-03, 1.484268e-03 }, 
    { 4.0265e+07, 376.15, 1.0, Fluid::N2,  1.909e-03, 1.851521e-03 }, 
    { 5.0209e+07, 376.15, 1.0, Fluid::N2,  2.236e-03, 2.180354e-03 }, 
    { 6.0348e+07, 376.15, 1.0, Fluid::N2,  2.585e-03, 2.483202e-03 }, 
    { 1.9986e+07, 376.15, 4.0, Fluid::N2,  5.128e-04, 5.270397e-04 }, 
    { 3.0125e+07, 376.15, 4.0, Fluid::N2,  7.194e-04, 7.379891e-04 }, 
    { 4.0070e+07, 376.15, 4.0, Fluid::N2,  8.832e-04, 9.195686e-04 }, 
    { 5.0014e+07, 376.15, 4.0, Fluid::N2,  1.033e-03, 1.081330e-03 }, 
    { 6.0153e+07, 376.15, 4.0, Fluid::N2,  1.189e-03, 1.229619e-03 }, 
    // H2S data from Selleck et al. (1952), Gillespie et al. (1984), Lee & Mather (1977)
    // Digitized from Fig 9 of Soreide & Whitson (1992)
    { 6.6480e+05, 311.15, 0.0, Fluid::H2S, 8.091e-03, 7.576558e-03 }, 
    { 1.0168e+06, 311.15, 0.0, Fluid::H2S, 1.219e-02, 1.176564e-02 }, 
    { 1.3687e+06, 311.15, 0.0, Fluid::H2S, 1.641e-02, 1.605412e-02 }, 
    { 1.7207e+06, 311.15, 0.0, Fluid::H2S, 2.063e-02, 2.044639e-02 }, 
    { 2.0726e+06, 311.15, 0.0, Fluid::H2S, 2.484e-02, 2.494013e-02 }, 
    { 2.6983e+06, 311.15, 0.0, Fluid::H2S, 3.316e-02, 3.308921e-02 }, 
    { 6.6480e+05, 344.15, 0.0, Fluid::H2S, 4.900e-03, 5.065083e-03 }, 
    { 1.0168e+06, 344.15, 0.0, Fluid::H2S, 7.407e-03, 7.913725e-03 }, 
    { 1.3687e+06, 344.15, 0.0, Fluid::H2S, 1.003e-02, 1.078026e-02 }, 
    { 1.7207e+06, 344.15, 0.0, Fluid::H2S, 1.265e-02, 1.366387e-02 }, 
    { 2.0726e+06, 344.15, 0.0, Fluid::H2S, 1.527e-02, 1.655986e-02 }, 
    { 2.7765e+06, 344.15, 0.0, Fluid::H2S, 2.040e-02, 2.237585e-02 }, 
    { 3.4609e+06, 344.15, 0.0, Fluid::H2S, 2.564e-02, 2.802407e-02 }, 
    { 3.9302e+06, 344.15, 0.0, Fluid::H2S, 3.054e-02, 3.186187e-02 }, 
    { 4.1453e+06, 344.15, 0.0, Fluid::H2S, 3.077e-02, 3.360257e-02 }, 
    { 4.2235e+06, 344.15, 0.0, Fluid::H2S, 3.248e-02, 3.423159e-02 }, 
    { 4.7709e+06, 344.15, 0.0, Fluid::H2S, 3.567e-02, 3.855773e-02 }, 
    { 4.8296e+06, 344.15, 0.0, Fluid::H2S, 3.624e-02, 3.901164e-02 }, 
    { 5.0642e+06, 344.15, 0.0, Fluid::H2S, 3.704e-02, 4.080062e-02 }, 
    { 5.2598e+06, 344.15, 0.0, Fluid::H2S, 3.932e-02, 4.225638e-02 }, 
    { 9.2877e+06, 344.15, 0.0, Fluid::H2S, 3.954e-02, 4.339253e-02 }, 
    { 1.0344e+07, 344.15, 0.0, Fluid::H2S, 3.977e-02, 4.374364e-02 }, 
    { 1.3804e+07, 344.15, 0.0, Fluid::H2S, 4.148e-02, 4.482798e-02 }, 
    { 1.3101e+06, 444.15, 0.0, Fluid::H2S, 2.849e-03, 2.370377e-03 }, 
    { 2.7179e+06, 444.15, 0.0, Fluid::H2S, 9.231e-03, 9.099828e-03 }, 
    { 4.1061e+06, 444.15, 0.0, Fluid::H2S, 1.538e-02, 1.574714e-02 }, 
    { 5.4944e+06, 444.15, 0.0, Fluid::H2S, 2.131e-02, 2.237576e-02 }, 
    { 6.8631e+06, 444.15, 0.0, Fluid::H2S, 2.724e-02, 2.885365e-02 }, 
    { 8.5838e+06, 444.15, 0.0, Fluid::H2S, 3.499e-02, 3.684163e-02 }, 
    { 1.0304e+07, 444.15, 0.0, Fluid::H2S, 4.342e-02, 4.453882e-02 }, 
    { 1.2045e+07, 444.15, 0.0, Fluid::H2S, 5.197e-02, 5.187989e-02 }, 
    // H2 data from Chabab et al. (2024)
    { 1.0001e+07, 298.15, 0.0, Fluid::H2,  1.360e-03, 1.359285e-03 }, 
    { 1.5001e+07, 298.15, 0.0, Fluid::H2,  1.997e-03, 1.995890e-03 }, 
    { 2.0001e+07, 298.15, 0.0, Fluid::H2,  2.639e-03, 2.608435e-03 }, 
    { 1.0111e+07, 323.40, 0.0, Fluid::H2,  1.251e-03, 1.270241e-03 }, 
    { 1.0131e+07, 323.40, 0.0, Fluid::H2,  1.239e-03, 1.272648e-03 }, 
    { 1.3001e+07, 323.40, 0.0, Fluid::H2,  1.593e-03, 1.614068e-03 }, 
    { 1.6501e+07, 323.40, 0.0, Fluid::H2,  2.011e-03, 2.020186e-03 }, 
    { 1.9991e+07, 323.40, 0.0, Fluid::H2,  2.442e-03, 2.414698e-03 }, 
    { 2.0011e+07, 323.40, 0.0, Fluid::H2,  2.455e-03, 2.416930e-03 }, 
    { 1.0001e+07, 373.15, 0.0, Fluid::H2,  1.412e-03, 1.380119e-03 }, 
    { 1.7511e+07, 373.15, 0.0, Fluid::H2,  2.433e-03, 2.355052e-03 }, 
    { 1.0071e+07, 298.15, 1.0, Fluid::H2,  1.070e-03, 1.136873e-03 }, 
    { 1.5001e+07, 298.15, 1.0, Fluid::H2,  1.603e-03, 1.657204e-03 }, 
    { 2.0001e+07, 298.15, 1.0, Fluid::H2,  2.136e-03, 2.164507e-03 }, 
    { 1.0031e+07, 323.40, 1.0, Fluid::H2,  1.021e-03, 1.066517e-03 }, 
    { 1.0061e+07, 323.40, 1.0, Fluid::H2,  1.021e-03, 1.069570e-03 }, 
    { 1.0101e+07, 323.40, 1.0, Fluid::H2,  1.024e-03, 1.073639e-03 }, 
    { 1.5006e+07, 323.40, 1.0, Fluid::H2,  1.514e-03, 1.562622e-03 }, 
    { 1.7501e+07, 323.40, 1.0, Fluid::H2,  1.767e-03, 1.804138e-03 }, 
    { 1.9991e+07, 323.40, 1.0, Fluid::H2,  2.035e-03, 2.040653e-03 }, 
    { 2.0001e+07, 323.40, 1.0, Fluid::H2,  2.035e-03, 2.041594e-03 }, 
    { 1.0011e+07, 373.15, 1.0, Fluid::H2,  1.197e-03, 1.201326e-03 }, 
    { 1.2601e+07, 373.15, 1.0, Fluid::H2,  1.485e-03, 1.498809e-03 }, 
    { 1.5036e+07, 373.15, 1.0, Fluid::H2,  1.776e-03, 1.772853e-03 }, 
    { 1.5046e+07, 373.15, 1.0, Fluid::H2,  1.776e-03, 1.773968e-03 }, 
    { 1.5071e+07, 373.15, 1.0, Fluid::H2,  1.773e-03, 1.776754e-03 }, 
    { 1.7551e+07, 373.15, 1.0, Fluid::H2,  2.040e-03, 2.050413e-03 }, 
    { 2.0046e+07, 373.15, 1.0, Fluid::H2,  2.346e-03, 2.320484e-03 }, 
    { 1.0001e+07, 298.15, 2.0, Fluid::H2,  8.864e-04, 9.339928e-04 }, 
    { 1.5001e+07, 298.15, 2.0, Fluid::H2,  1.322e-03, 1.369797e-03 }, 
    { 2.0001e+07, 298.15, 2.0, Fluid::H2,  1.718e-03, 1.788160e-03 }, 
    { 1.0001e+07, 323.40, 2.0, Fluid::H2,  8.826e-04, 8.954498e-04 }, 
    { 1.5001e+07, 323.40, 2.0, Fluid::H2,  1.301e-03, 1.314657e-03 }, 
    { 2.0001e+07, 323.40, 2.0, Fluid::H2,  1.724e-03, 1.717318e-03 }, 
    { 1.0001e+07, 373.15, 2.0, Fluid::H2,  9.938e-04, 1.038365e-03 }, 
    { 1.5001e+07, 373.15, 2.0, Fluid::H2,  1.519e-03, 1.529530e-03 }, 
    { 2.0051e+07, 373.15, 2.0, Fluid::H2,  2.050e-03, 2.005806e-03 }, 
    { 1.0001e+07, 298.15, 4.0, Fluid::H2,  5.942e-04, 6.345449e-04 }, 
    { 1.5001e+07, 298.15, 4.0, Fluid::H2,  9.359e-04, 9.297604e-04 }, 
    { 2.0001e+07, 298.15, 4.0, Fluid::H2,  1.218e-03, 1.212645e-03 }, 
    { 1.0001e+07, 323.40, 4.0, Fluid::H2,  6.174e-04, 6.304264e-04 }, 
    { 1.5001e+07, 323.40, 4.0, Fluid::H2,  9.575e-04, 9.247533e-04 }, 
    { 2.0001e+07, 323.40, 4.0, Fluid::H2,  1.292e-03, 1.206997e-03 }, 
    { 1.0001e+07, 373.15, 4.0, Fluid::H2,  7.799e-04, 7.714822e-04 }, 
    { 1.5001e+07, 373.15, 4.0, Fluid::H2,  1.145e-03, 1.135169e-03 }, 
    { 2.0001e+07, 373.15, 4.0, Fluid::H2,  1.575e-03, 1.483850e-03 }, 
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
