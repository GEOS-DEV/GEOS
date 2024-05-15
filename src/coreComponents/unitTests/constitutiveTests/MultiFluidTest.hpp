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
 * @file MultiFluidTest.hpp
 */

#ifndef GEOS_UNITTESTS_CONSTITUTIVETESTS_MULTIFLUIDTEST_HPP_
#define GEOS_UNITTESTS_CONSTITUTIVETESTS_MULTIFLUIDTEST_HPP_

#include "constitutiveTestHelpers.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

using namespace geos::dataRepository;

namespace geos
{
namespace testing
{

template< typename FLUID_TYPE >
struct UsePVTPackage
{
  static constexpr bool value = false;
};

template< integer N >
using Feed = std::array< real64, N >;

template< integer N, typename ARRAY1, typename ARRAY2 >
void copy( ARRAY1 const & inputArray, ARRAY2 & outputArray )
{
  for( integer i = 0; i < N; ++i )
  {
    outputArray[i] = inputArray[i];
  }
}

template< integer N, typename T >
void fill( array1d< T > & outputArray, std::array< T const, N > const inputArray )
{
  for( integer i = 0; i < N; ++i )
  {
    outputArray.emplace_back( inputArray[i] );
  }
}

template< integer N, integer U, typename T >
void fill( arraySlice1d< T, U > const & outputArray, std::array< T const, N > const inputArray )
{
  for( integer i = 0; i < N; ++i )
  {
    outputArray[i] = inputArray[i];
  }
}

template< integer NUM_PHASE, integer NUM_COMP >
struct MultiFluidTestData
{
  MultiFluidTestData( real64 const pressure_,
                      real64 const temperature_,
                      Feed< NUM_COMP > && composition_ );

  MultiFluidTestData( real64 const pressure_,
                      real64 const temperature_,
                      arraySlice1d< real64 const > const & composition_ );

  real64 const pressure{0.0};
  real64 const temperature{0.0};
  Feed< NUM_COMP > composition{};
  real64 const totalDensity{0.0};
  Feed< NUM_PHASE > phaseFraction{};
  Feed< NUM_PHASE > phaseDensity{};
  Feed< NUM_PHASE > phaseMassDensity{};
  Feed< NUM_PHASE > phaseViscosity{};
};

template< typename FLUID_TYPE, integer NUM_PHASE, integer NUM_COMP >
class MultiFluidTest : public ConstitutiveTestBase< MultiFluidBase >
{
public:
  static constexpr integer numPhase = NUM_PHASE;
  static constexpr integer numComp = NUM_COMP;
  using TestData = MultiFluidTestData< numPhase, numComp >;
  using PhaseFeed = Feed< numPhase >;
  using CompositionFeed = Feed< numComp >;

public:
  MultiFluidTest() = default;
  ~MultiFluidTest() override = default;

  MultiFluidBase & getFluid() const { return *m_model; }

  Group & getParent() { return m_parent; }

  void testValuesAgainstPreviousImplementation( typename FLUID_TYPE::KernelWrapper const & wrapper,
                                                MultiFluidTestData< NUM_PHASE, NUM_COMP > const & testData,
                                                real64 const relTol ) const;

  void testNumericalDerivatives( MultiFluidBase & fluid,
                                 Group * parent,
                                 MultiFluidTestData< NUM_PHASE, NUM_COMP > const & testData,
                                 real64 const perturbParameter,
                                 real64 const relTol,
                                 real64 const absTol = std::numeric_limits< real64 >::max() );


protected:
  virtual void resetFluid( MultiFluidBase & fluid ) const
  {
    GEOS_UNUSED_VAR( fluid );
  }

  static void writeTableToFile( string const & fileName, char const * str )
  {
    std::ofstream os( fileName );
    ASSERT_TRUE( os.is_open() );
    os << str;
    os.close();
  }

  static void removeFile( string const & fileName )
  {
    int const ret = std::remove( fileName.c_str() );
    ASSERT_TRUE( ret == 0 );
  }
};

template< typename FLUID_TYPE, integer NUM_PHASE, integer NUM_COMP >
void MultiFluidTest< FLUID_TYPE, NUM_PHASE, NUM_COMP >::
testNumericalDerivatives( MultiFluidBase & fluid,
                          Group * parent,
                          MultiFluidTestData< NUM_PHASE, NUM_COMP > const & testData,
                          real64 const perturbParameter,
                          real64 const relTol,
                          real64 const absTol )
{
  using Deriv = multifluid::DerivativeOffset;

  integer const NC = fluid.numFluidComponents();
  integer const NP = fluid.numFluidPhases();
  integer const NDOF = NC+2;

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, NC );
  for( integer i = 0; i < NC; ++i )
  {
    compositionValues[0][i] = testData.composition[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  auto const & components = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  auto const & phases     = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );

  // create a clone of the fluid to run updates on
  std::unique_ptr< ConstitutiveBase > fluidCopyPtr = fluid.deliverClone( "fluidCopy", parent );
  MultiFluidBase & fluidCopy = dynamicCast< MultiFluidBase & >( *fluidCopyPtr );

  fluid.allocateConstitutiveData( fluid.getParent(), 1 );
  fluidCopy.allocateConstitutiveData( fluid.getParent(), 1 );

  // extract data views from both fluids
  #define GET_FLUID_DATA( FLUID, TRAIT ) \
    FLUID.getReference< TRAIT::type >( TRAIT::key() )[0][0]

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseFrac {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseFraction ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseFraction )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseDens {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseDensity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseDensity )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseVisc {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseViscosity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseViscosity )
  };

  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 > phaseCompFrac {
    GET_FLUID_DATA( fluid, fields::multifluid::phaseCompFraction ),
    GET_FLUID_DATA( fluid, fields::multifluid::dPhaseCompFraction )
  };

  MultiFluidVarSlice< real64, 0, USD_FLUID - 2, USD_FLUID_DC - 2 > totalDens {
    GET_FLUID_DATA( fluid, fields::multifluid::totalDensity ),
    GET_FLUID_DATA( fluid, fields::multifluid::dTotalDensity )
  };

  auto const & phaseFracCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseFraction );
  auto const & phaseDensCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseDensity );
  auto const & phaseViscCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseViscosity );
  auto const & phaseCompFracCopy = GET_FLUID_DATA( fluidCopy, fields::multifluid::phaseCompFraction );
  auto const & totalDensCopy     = GET_FLUID_DATA( fluidCopy, fields::multifluid::totalDensity );

#undef GET_FLUID_DATA

  real64 const pressure = testData.pressure;
  real64 const temperature = testData.temperature;

  // set the original fluid state to current
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    fluidWrapper.update( 0, 0, pressure, temperature, composition );
  } );

  // now perturb variables and update the copied fluid's state
  constitutive::constitutiveUpdatePassThru( fluidCopy, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // to be able to use the checkDerivative utility function, we have to invert the layout
    auto dPhaseFrac     = invertLayout( phaseFrac.derivs.toSliceConst(), NP, NDOF );
    auto dPhaseDens     = invertLayout( phaseDens.derivs.toSliceConst(), NP, NDOF );
    auto dPhaseVisc     = invertLayout( phaseVisc.derivs.toSliceConst(), NP, NDOF );
    auto dTotalDens     = invertLayout( totalDens.derivs.toSliceConst(), NDOF );
    auto dPhaseCompFrac = invertLayout( phaseCompFrac.derivs.toSliceConst(), NP, NC, NDOF );

    // update pressure and check derivatives
    {
      real64 const dP = perturbParameter * (pressure + perturbParameter);
      resetFluid( fluidCopy );
      fluidWrapper.update( 0, 0, pressure + dP, temperature, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseFrac", "Pres", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseDens", "Pres", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseVisc", "Pres", phases );
      checkDerivative( totalDensCopy, totalDens.value, dTotalDens[Deriv::dP],
                       dP, relTol, absTol, "totalDens", "Pres" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac[Deriv::dP].toSliceConst(),
                       dP, relTol, absTol, "phaseCompFrac", "Pres", phases, components );
    }

    // update temperature and check derivatives
    {
      real64 const dT = perturbParameter * (temperature + perturbParameter);
      resetFluid( fluidCopy );
      fluidWrapper.update( 0, 0, pressure, temperature + dT, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseFrac", "Temp", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseDens", "Temp", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseVisc", "Temp", phases );
      checkDerivative( totalDensCopy, totalDens.value, dTotalDens[Deriv::dT],
                       dT, relTol, absTol, "totalDens", "Temp" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac[Deriv::dT].toSliceConst(),
                       dT, relTol, absTol, "phaseCompFrac", "Temp", phases, components );
    }

    array2d< real64, compflow::LAYOUT_COMP > compNew( 1, NC );
    for( integer jc = 0; jc < NC; ++jc )
    {
      real64 const dC = perturbParameter * ( composition[jc] + perturbParameter );
      for( integer ic = 0; ic < NC; ++ic )
      {
        compNew[0][ic] = composition[ic];
      }
      compNew[0][jc] += dC;

      // Note: in PVTPackage, derivatives are obtained with finite-difference approx **with normalization of the comp fraction**
      //       The component fraction is perturbed (just as above), and then all the component fractions are normalized (as below)
      //       But, in the native DO model and in CO2BrinePhillips, derivatives are computed analytically, which results in different
      //       derivatives wrt component fractions--although the derivatives wrt component densities obtained with the chain rule
      //       in the solver will be very similar (see discussion on PR #1325 on GitHub).
      //
      //       Since both approaches--FD approximation of derivatives with normalization, and analytical derivatives--are correct,
      //       we have to support both when we check the intermediate derivatives wrt component fractions below. Therefore, if the
      //       PVTPackage is used, then we normalize the perturbed component fractions before taking the FD approx. If the native
      //       DO or CO2-brine models are used, we skip the normalization below.
      if constexpr ( UsePVTPackage< FLUID_TYPE >::value )
      {
        // renormalize
        real64 sum = 0.0;
        for( integer ic = 0; ic < NC; ++ic )
        {
          sum += compNew[0][ic];
        }
        for( integer ic = 0; ic < NC; ++ic )
        {
          compNew[0][ic] /= sum;
        }
      }

      resetFluid( fluidCopy );
      fluidWrapper.update( 0, 0, pressure, temperature, compNew[0] );

      string const var = "compFrac[" + components[jc] + "]";
      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac[Deriv::dC+jc].toSliceConst(),
                       dC, relTol, absTol, "phaseFrac", var, phases );/**
                                                                         checkDerivative( phaseDensCopy.toSliceConst(),
                                                                            phaseDens.value.toSliceConst(),
                                                                            dPhaseDens[Deriv::dC+jc].toSliceConst(),
                                                                         dC, relTol, absTol, "phaseDens", var, phases );
                                                                         checkDerivative( phaseViscCopy.toSliceConst(),
                                                                            phaseVisc.value.toSliceConst(),
                                                                            dPhaseVisc[Deriv::dC+jc].toSliceConst(),
                                                                         dC, relTol, absTol, "phaseVisc", var, phases );
                                                                         checkDerivative( totalDensCopy, totalDens.value,
                                                                            dTotalDens[Deriv::dC+jc],
                                                                         dC, relTol, absTol, "totalDens", var );
                                                                         checkDerivative( phaseCompFracCopy.toSliceConst(),
                                                                            phaseCompFrac.value.toSliceConst(),
                                                                            dPhaseCompFrac[Deriv::dC+jc].toSliceConst(),
                                                                         dC, relTol, absTol, "phaseCompFrac", var, phases, components );**/
    }
  } );
}

template< typename FLUID_TYPE, integer NUM_PHASE, integer NUM_COMP >
void MultiFluidTest< FLUID_TYPE, NUM_PHASE, NUM_COMP >::testValuesAgainstPreviousImplementation( typename FLUID_TYPE::KernelWrapper const & wrapper,
                                                                                                 MultiFluidTestData< NUM_PHASE, NUM_COMP > const & testData,
                                                                                                 real64 const relTol ) const
{
  integer constexpr numDof = numComp + 2;

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, numComp );
  for( integer i = 0; i < numComp; ++i )
  {
    compositionValues[0][i] = testData.composition[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFraction( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFraction( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseDensity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseDensity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseMassDensity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseMassDensity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseViscosity( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseViscosity( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseEnthalpy( 1, 1, numPhase, numDof );
  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseInternalEnergy( 1, 1, numPhase, numDof );
  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFraction( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFraction( 1, 1, numPhase, numComp, numDof );
  StackArray< real64, 2, 1, LAYOUT_FLUID > totalDensity( 1, 1 );
  StackArray< real64, 3, numDof, LAYOUT_FLUID_DC >  dTotalDensity( 1, 1, numDof );

  wrapper.compute( testData.pressure,
                   testData.temperature,
                   composition,
                   { phaseFraction[0][0], dPhaseFraction[0][0] },
                   { phaseDensity[0][0], dPhaseDensity[0][0] },
                   { phaseMassDensity[0][0], dPhaseMassDensity[0][0] },
                   { phaseViscosity[0][0], dPhaseViscosity[0][0] },
                   { phaseEnthalpy[0][0], dPhaseEnthalpy[0][0] },
                   { phaseInternalEnergy[0][0], dPhaseInternalEnergy[0][0] },
                   { phaseCompFraction[0][0], dPhaseCompFraction[0][0] },
                   { totalDensity[0][0], dTotalDensity[0][0] } );

  checkRelativeError( totalDensity[0][0], testData.totalDensity, relTol );
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    checkRelativeError( phaseFraction[0][0][ip], testData.phaseFraction[ip], relTol );
    checkRelativeError( phaseDensity[0][0][ip], testData.phaseDensity[ip], relTol );
    checkRelativeError( phaseMassDensity[0][0][ip], testData.phaseMassDensity[ip], relTol );
    checkRelativeError( phaseViscosity[0][0][ip], testData.phaseViscosity[ip], relTol );
  }
}

template< integer NUM_PHASE, integer NUM_COMP >
MultiFluidTestData< NUM_PHASE, NUM_COMP >::MultiFluidTestData( real64 const pressure_,
                                                               real64 const temperature_,
                                                               Feed< NUM_COMP > && composition_ ):
  pressure( pressure_ ),
  temperature( temperature_ ),
  composition( composition_ )
{}

template< integer NUM_PHASE, integer NUM_COMP >
MultiFluidTestData< NUM_PHASE, NUM_COMP >::MultiFluidTestData( real64 const pressure_,
                                                               real64 const temperature_,
                                                               arraySlice1d< real64 const > const & composition_ ):
  pressure( pressure_ ),
  temperature( temperature_ )
{
  copy< NUM_COMP >( composition_, composition );
}

} // namespace testing

} // namespace geos

#endif //GEOS_UNITTESTS_CONSTITUTIVETESTS_MULTIFLUIDTEST_HPP_
