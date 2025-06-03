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
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp"
#include "TestFluid.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

/**
 * This test compares the solubility calculated from our implementation of the Soreide Whitson
 * model against published experimental data.
 */
template< int NCOMP >
using SoreideWhitsonFlashData = std::tuple<
  real64 const,     // pressure
  real64 const,     // temperature
  Feed< NCOMP > const,// composition
  real64 const      // liquid mole fraction
  >;

template< int NCOMP >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    auto fluid = TestFluid< 4 >::create( {0, 0, 0, 0} );
    fluid->componentNames = { "H2", "CH4", "CO2", "H2O" };
    TestFluid< 4 >::populateArray( fluid->criticalPressure, Feed< 4 >{1.29640e+06, 4.59920e+06, 7.37730e+06, 2.20640e+07} );
    TestFluid< 4 >::populateArray( fluid->criticalTemperature, Feed< 4 >{3.31450e+01, 1.90564e+02, 3.04128e+02, 6.47096e+02} );
    TestFluid< 4 >::populateArray( fluid->criticalVolume, Feed< 4 >{6.44828e-05, 9.86278e-05, 9.41185e-05, 5.59480e-05} );
    TestFluid< 4 >::populateArray( fluid->acentricFactor, Feed< 4 >{-2.19000e-01, 1.14200e-02, 2.23940e-01, 3.44300e-01} );
    TestFluid< 4 >::populateArray( fluid->molecularWeight, Feed< 4 >{2.01588e-03, 1.60425e-02, 4.40095e-02, 1.80153e-02} );
    fluid->setBinaryCoefficients( Feed< 6 >{ 0.0000, 0.0000, 0.0000, -0.3776, 0.4850, 0.1896 } );
    return fluid;
  }
};

template< integer NCOMP, integer PPM=0 >
class SoreideWhitsonFlashMultiComponentTestFixture : public ::testing::TestWithParam< SoreideWhitsonFlashData< NCOMP > >
{
protected:
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numPhases = 2;
  static constexpr int numComps = NCOMP;
  static constexpr int numDofs = NCOMP + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using PhasePropSlice = NegativeTwoPhaseFlashModelUpdate::PhaseProp::SliceType;
  using PhaseCompSlice = NegativeTwoPhaseFlashModelUpdate::PhaseComp::SliceType;

public:
  SoreideWhitsonFlashMultiComponentTestFixture():
    m_fluid( FluidData< NCOMP >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = NegativeTwoPhaseFlashModel::createParameters( std::move( m_parameters ) );
    m_parameters = BrineSalinity::create( std::move( m_parameters ) );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    equationOfState->m_equationsOfStateNames = {
      EnumStrings< EquationOfStateType >::toString( EquationOfStateType::SoreideWhitson ),
      EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson )
    };

    auto * brineSalinity = const_cast< BrineSalinity * >(m_parameters->get< BrineSalinity >());
    real64 const massFraction = 1.0e-6*PPM;
    brineSalinity->m_salinity = massFraction / 58.44e-3;

    m_flash = std::make_unique< NegativeTwoPhaseFlashModel >( "Flash", componentProperties, *m_parameters );
  }

  ~SoreideWhitsonFlashMultiComponentTestFixture() = default;

  void testFlash( SoreideWhitsonFlashData< NCOMP > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NCOMP >::createArray( composition, std::get< 2 >( data ));

    stackArray2d< real64, numComps > kValues( 1, numComps );
    kValues( 0, 0 ) = 0.0;

    StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
    StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
    StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
    StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

    auto phaseFraction = phaseFractionData[0][0];
    auto dPhaseFraction = dPhaseFractionData[0][0];
    auto phaseComponentFraction = phaseComponentFractionData[0][0];
    auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

    auto componentProperties = m_fluid->createKernelWrapper();
    auto flashKernelWrapper = m_flash->createKernelWrapper();

    flashKernelWrapper.compute( componentProperties,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                kValues.toSlice(),
                                PhasePropSlice( phaseFraction, dPhaseFraction ),
                                PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

std::cout << phaseFraction.toSliceConst() << std::endl;

  }

protected:
  std::unique_ptr< TestFluid< NCOMP > > m_fluid{};
  std::unique_ptr< NegativeTwoPhaseFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using SoreideWhitson4 = SoreideWhitsonFlashMultiComponentTestFixture< 4 >;

TEST_P( SoreideWhitson4, testFlash )
{
  testFlash( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonFlashMultiComponent, SoreideWhitson4,
  ::testing::ValuesIn< SoreideWhitsonFlashData<4> >({
    { 1.00000e+05, 333.15, {0.4, 0.1, 0.3, 0.2}, 0.0 }
  })
);

/* UNCRUSTIFY-ON */

} // testing
} // geos
