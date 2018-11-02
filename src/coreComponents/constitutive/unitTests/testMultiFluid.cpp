/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include <constitutive/Fluid/CompositionalMultiphaseFluid.hpp>
#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/Fluid/CompositionalMultiphaseFluid.hpp"
#include "constitutive/Fluid/BlackOilFluid.hpp"

using namespace geosx;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

template<typename T, int NDIM>
using array = multidimensionalArray::ManagedArray<T,NDIM,localIndex>;

// helper struct to represent a var and its derivatives (always with array views, not pointers)
template<int DIM>
struct TestCompositionalVarContainer
{
  array_view<real64,DIM>   value; // variable value
  array_view<real64,DIM>   dPres; // derivative w.r.t. pressure
  array_view<real64,DIM>   dTemp; // derivative w.r.t. temperature
  array_view<real64,DIM+1> dComp; // derivative w.r.t. composition
};

inline bool checkRelativeError( real64 v1, real64 v2, real64 relTol, string const & name )
{
  real64 const delta = std::fabs( v1 - v2 );
  real64 const value = std::fmax( std::fabs(v1), std::fabs(v2) );
  if (value > 0)
  {
    real64 const err = delta / value;
    if (err > relTol)
    {
      GEOS_LOG("[ " << name << " ] relative error: " << err << "(" << v1 << " vs " << v2 << ")");
      return false;
    }
  }
  return true;
}

bool checkDerivative( real64 const & valueEps, real64 const & value, real64 const & deriv,
                      real64 eps, real64 relTol,
                      string const & name, string const & var )
{
  real64 const numDeriv = (valueEps - value) / eps;
  return checkRelativeError( numDeriv, deriv, relTol, "d(" + name + ")/d(" + var + ")" );
}

template<int DIM, typename ... Args>
typename std::enable_if<DIM != 0, bool>::type
checkDerivative( array_view<real64,DIM> const & valueEps,
                 array_view<real64,DIM> const & value,
                 array_view<real64,DIM> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists)
{
  const auto size = valueEps.size(0);

  bool result = true;
  for (localIndex i = 0; i < size; ++i)
  {
    result &= checkDerivative( valueEps.slice(i), value.slice(i), deriv.slice(i), eps, relTol,
                               name + "[" + labels[i] + "]", var, label_lists... );
  }

  return result;
}

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
template<int DIM>
array<real64,DIM> invertLayout( array_view<real64,DIM> const & input )
{
  array<real64,DIM> output(input);
  return output;
}

template<>
array<real64,2> invertLayout( array_view<real64,2> const & input )
{
  array<real64,2> output(input.size(1), input.size(0));

  for (int i = 0; i < input.size(0); ++i)
    for (int j = 0; j < input.size(1); ++j)
      output[j][i] = input[i][j];

  return output;
}

template<>
array<real64,3> invertLayout( array_view<real64,3> const & input )
{
  array<real64,3> output(input.size(2), input.size(0), input.size(1));

  for (int i = 0; i < input.size(0); ++i)
    for (int j = 0; j < input.size(1); ++j)
      for (int k = 0; k < input.size(2); ++k)
      output[k][i][j] = input[i][j][k];

  return output;
}

bool testNumericalDerivatives( MultiFluidBase * fluid,
                               real64 P,
                               real64 T,
                               array1d<real64> const & composition,
                               real64 perturbParameter = 1e-6,
                               real64 relTol = 1e-4 )
{
  localIndex const NC = fluid->numFluidComponents();

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  // create a clone of the fluid to run updates on
  auto fluidCopy = fluid->DeliverClone( "fluidCopy", nullptr );

  fluid->AllocateConstitutiveData( fluid->getParent(), 1 );
  fluidCopy->AllocateConstitutiveData( fluid->getParent(), 1 );

  // extract data views from both fluids
#define GET_FLUID_DATA( FLUID, DIM, KEY ) \
  FLUID->getReference<array<real64,DIM>>( MultiFluidBase::viewKeyStruct::KEY ).slice(0,0)

  TestCompositionalVarContainer<1> phaseFrac {
    GET_FLUID_DATA( fluid, 3, phaseFractionString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseFraction_dGlobalCompFractionString )
  };

  TestCompositionalVarContainer<1> phaseDens {
    GET_FLUID_DATA( fluid, 3, phaseDensityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseDensity_dGlobalCompFractionString )
  };

  TestCompositionalVarContainer<1> phaseVisc {
    GET_FLUID_DATA( fluid, 3, phaseViscosityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseViscosity_dGlobalCompFractionString )
  };

  TestCompositionalVarContainer<2> phaseCompFrac {
    GET_FLUID_DATA( fluid, 4, phaseCompFractionString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 5, dPhaseCompFraction_dGlobalCompFractionString )
  };

  TestCompositionalVarContainer<0> totalDens {
    GET_FLUID_DATA( fluid, 2, totalDensityString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 3, dTotalDensity_dGlobalCompFractionString )
  };

  auto const & phaseFracCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseFractionString );
  auto const & phaseDensCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseDensityString );
  auto const & phaseViscCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseViscosityString );
  auto const & phaseCompFracCopy = GET_FLUID_DATA( fluidCopy, 4, phaseCompFractionString );
  auto const & totalDensCopy     = GET_FLUID_DATA( fluidCopy, 2, totalDensityString );

#undef GET_FLUID_DATA

  // extract data views (without derivatives) from the copy

  // set the fluid state to current
  fluid->StateUpdatePointMultiphaseFluid( P, T, composition.data(), 0, 0 );

  bool result = true;

  // update pressure and check derivatives
  {
    real64 const dP = perturbParameter * (P + perturbParameter);
    fluidCopy->StateUpdatePointMultiphaseFluid( P + dP, T, composition.data(), 0, 0 );

    result &= checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dPres, dP, relTol, "phaseFrac", "P", phases );
    result &= checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dPres, dP, relTol, "phaseDens", "P", phases );
    result &= checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dPres, dP, relTol, "phaseVisc", "P", phases );
    result &= checkDerivative( totalDensCopy, totalDens.value, totalDens.dPres, dP, relTol, "totalDens", "P" );
    result &= checkDerivative( phaseCompFracCopy, phaseCompFrac.value, phaseCompFrac.dPres,
                               dP, relTol, "phaseCompFrac", "P", phases, components );
  }

  // update temperature and check derivatives
  {
    real64 const dT = perturbParameter * (T + perturbParameter);
    fluidCopy->StateUpdatePointMultiphaseFluid( P, T + dT, composition.data(), 0, 0 );

    result &= checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dTemp, dT, relTol, "phaseFrac", "T", phases );
    result &= checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dTemp, dT, relTol, "phaseDens", "T", phases );
    result &= checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dTemp, dT, relTol, "phaseVisc", "T", phases );
    result &= checkDerivative( totalDensCopy, totalDens.value, totalDens.dTemp, dT, relTol, "totalDens", "T" );
    result &= checkDerivative( phaseCompFracCopy, phaseCompFrac.value, phaseCompFrac.dTemp,
                               dT, relTol, "phaseCompFrac", "T", phases, components );
  }

  // update temperature and check derivatives
  auto dPhaseFrac_dC     = invertLayout(phaseFrac.dComp);
  auto dPhaseDens_dC     = invertLayout(phaseDens.dComp);
  auto dPhaseVisc_dC     = invertLayout(phaseVisc.dComp);
  auto dTotalDens_dC     = invertLayout(totalDens.dComp);
  auto dPhaseCompFrac_dC = invertLayout(phaseCompFrac.dComp);

  array1d<real64> compNew;
  for (localIndex jc = 0; jc < NC; ++jc)
  {
    real64 const dC = perturbParameter * (composition[jc] + perturbParameter);
    compNew = composition;
    compNew[jc] += dC;

    // renormalize
    real64 sum = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
      sum += compNew[ic];
    for (localIndex ic = 0; ic < NC; ++ic)
      compNew[ic] /= sum;

    fluidCopy->StateUpdatePointMultiphaseFluid( P, T, compNew.data(), 0, 0 );
    string var = "C[" + components[jc] + "]";

    result &= checkDerivative( phaseFracCopy, phaseFrac.value, dPhaseFrac_dC.slice(jc), dC, relTol, "phaseFrac", var, phases );
    result &= checkDerivative( phaseDensCopy, phaseDens.value, dPhaseDens_dC.slice(jc), dC, relTol, "phaseDens", var, phases );
    result &= checkDerivative( phaseViscCopy, phaseVisc.value, dPhaseVisc_dC.slice(jc), dC, relTol, "phaseVisc", var, phases );
    result &= checkDerivative( totalDensCopy, totalDens.value, dTotalDens_dC.slice(jc), dC, relTol, "totalDens", var );
    result &= checkDerivative( phaseCompFracCopy, phaseCompFrac.value, dPhaseCompFrac_dC.slice(jc),
                               dC, relTol, "phaseCompFrac", var, phases, components );
  }

  return result;
}

MultiFluidBase * makeCompositionalFluid( string const & name, ManagedGroup * parent )
{
  auto fluid = parent->RegisterGroup<CompositionalMultiphaseFluid>( name );

  auto & phaseNames = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & eqnOfState = fluid->getReference<string_array>( CompositionalMultiphaseFluid::viewKeyStruct::equationsOfStateString );
  eqnOfState.resize( 2 );
  eqnOfState[0] = "PR"; eqnOfState[1] = "PR";

  auto & compNames = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  compNames.resize( 4 );
  compNames[0] = "N2"; compNames[1] = "C10"; compNames[2] = "C20"; compNames[3] = "H20";

  auto & critPres = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalPressureString );
  critPres.resize( 4 );
  critPres[0] = 34e5; critPres[1] = 25.3e5; critPres[2] = 14.6e5; critPres[3] = 220.5e5;

  auto & critTemp = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalTemperatureString );
  critTemp.resize( 4 );
  critTemp[0] = 126.2; critTemp[1] = 622.0; critTemp[2] = 782.0; critTemp[3] = 647.0;

  auto & acFactor = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentAcentricFactorString );
  acFactor.resize( 4 );
  acFactor[0] = 0.04; acFactor[1] = 0.443; acFactor[2] = 0.816; acFactor[3] = 0.344;

  auto & molarWgt = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentMolarWeightString );
  molarWgt.resize( 4 );
  molarWgt[0] = 28e-3; molarWgt[1] = 134e-3; molarWgt[2] = 275e-3; molarWgt[3] = 18e-3;

  return fluid;
}

TEST(testMultiFluid, numericalDerivativesCompositionalFluid)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  MultiFluidBase * fluid = makeCompositionalFluid( "fluid", parent.get() );

  parent->Initialize( parent.get() );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d<real64> comp(4);
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  auto sqrtprecision = sqrt(std::numeric_limits<real64>::epsilon());

  bool derivsAreOk = testNumericalDerivatives( fluid, P, T, comp, sqrtprecision, 1e-3 );

  EXPECT_TRUE( derivsAreOk );
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI

  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}
