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
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
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
using array = LvArray::Array<T,NDIM,localIndex>;

template<typename T>
::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *,
                                                     T v1, T v2, T relTol )
{
  T const delta = std::abs( v1 - v2 );
  T const value = std::max( std::abs(v1), std::abs(v2) );
  if (delta > relTol * value)
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision(5)
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
  }
  return ::testing::AssertionSuccess();
}

template<typename T>
void checkRelativeError( T v1, T v2, T relTol )
{
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

template<typename T>
void checkRelativeError( T v1, T v2, T relTol, string const & name )
{
  SCOPED_TRACE(name);
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

template<typename T>
void checkDerivative( T valueEps, T value, T deriv, real64 eps, real64 relTol, string const & name, string const & var )
{
  T numDeriv = (valueEps - value) / eps;
  checkRelativeError( deriv, numDeriv, relTol, "d(" + name + ")/d(" + var + ")" );
}

template<typename T, typename ... Args>
void
checkDerivative( arraySlice1d<T> const & valueEps,
                 arraySlice1d<T> const & value,
                 arraySlice1d<T> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  localIndex const size = labels.size(0);

  for (localIndex i = 0; i < size; ++i)
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

template<typename T, int DIM, typename ... Args>
typename std::enable_if<(DIM > 1), void>::type
checkDerivative( array_slice<T,DIM> const & valueEps,
                 array_slice<T,DIM> const & value,
                 array_slice<T,DIM> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  localIndex const size = labels.size(0);

  for (localIndex i = 0; i < size; ++i)
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
array1d<real64> invertLayout( arraySlice1d<real64 const> const & input, localIndex N )
{
  array<real64,1> output( N );
  for (int i = 0; i < N; ++i)
    output[i] = input[i];

  return output;
}

array2d<real64> invertLayout( arraySlice2d<real64 const> const & input, localIndex N1, localIndex N2 )
{
  array<real64,2> output( N2, N1 );

  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j)
      output[j][i] = input[i][j];

  return output;
}

array3d<real64> invertLayout( arraySlice3d<real64 const> const & input, localIndex N1, localIndex N2, localIndex N3 )
{
  array<real64,3> output( N3, N1, N2 );

  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j)
      for (int k = 0; k < N3; ++k)
        output[k][i][j] = input[i][j][k];

  return output;
}

void testNumericalDerivatives( MultiFluidBase * fluid,
                               real64 P,
                               real64 T,
                               arraySlice1d<real64> const & composition,
                               real64 perturbParameter,
                               real64 relTol )
{
  localIndex const NC = fluid->numFluidComponents();
  localIndex const NP = fluid->numFluidPhases();

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  // create a clone of the fluid to run updates on
  auto fluidCopyPtr = fluid->DeliverClone( "fluidCopy", nullptr );
  auto fluidCopy = fluidCopyPtr->group_cast<MultiFluidBase *>();

  fluid->AllocateConstitutiveData( fluid->getParent(), 1 );
  fluidCopy->AllocateConstitutiveData( fluid->getParent(), 1 );

  // extract data views from both fluids
#define GET_FLUID_DATA( FLUID, DIM, KEY ) \
  FLUID->getReference<array<real64,DIM>>( MultiFluidBase::viewKeyStruct::KEY )[0][0]

  CompositionalVarContainer<1> phaseFrac {
    GET_FLUID_DATA( fluid, 3, phaseFractionString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseFraction_dGlobalCompFractionString )
  };

  CompositionalVarContainer<1> phaseDens {
    GET_FLUID_DATA( fluid, 3, phaseDensityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseDensity_dGlobalCompFractionString )
  };

  CompositionalVarContainer<1> phaseVisc {
    GET_FLUID_DATA( fluid, 3, phaseViscosityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseViscosity_dGlobalCompFractionString )
  };

  CompositionalVarContainer<2> phaseCompFrac {
    GET_FLUID_DATA( fluid, 4, phaseCompFractionString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 5, dPhaseCompFraction_dGlobalCompFractionString )
  };

  CompositionalVarContainer<0> totalDens {
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

  // set the fluid state to current
  fluid->StateUpdatePointMultiFluid(P, T, composition, 0, 0);

  // update pressure and check derivatives
  {
    real64 const dP = perturbParameter * (P + perturbParameter);
    fluidCopy->StateUpdatePointMultiFluid( P + dP, T, composition, 0, 0 );

    checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dPres, dP, relTol, "phaseFrac", "Pres", phases );
    checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dPres, dP, relTol, "phaseDens", "Pres", phases );
    checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dPres, dP, relTol, "phaseVisc", "Pres", phases );
    checkDerivative( totalDensCopy, totalDens.value, totalDens.dPres, dP, relTol, "totalDens", "Pres" );
    checkDerivative( phaseCompFracCopy, phaseCompFrac.value, phaseCompFrac.dPres, dP, relTol,
                     "phaseCompFrac", "Pres", phases, components );
  }

  // update temperature and check derivatives
  {
    real64 const dT = perturbParameter * (T + perturbParameter);
    fluidCopy->StateUpdatePointMultiFluid( P, T + dT, composition, 0, 0 );

    checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dTemp, dT, relTol, "phaseFrac", "Temp", phases );
    checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dTemp, dT, relTol, "phaseDens", "Temp", phases );
    checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dTemp, dT, relTol, "phaseVisc", "Temp", phases );
    checkDerivative( totalDensCopy, totalDens.value, totalDens.dTemp, dT, relTol, "totalDens", "Temp" );
    checkDerivative( phaseCompFracCopy, phaseCompFrac.value, phaseCompFrac.dTemp, dT, relTol,
                     "phaseCompFrac", "Temp", phases, components );
  }

  // update composition and check derivatives
  auto dPhaseFrac_dC     = invertLayout( phaseFrac.dComp, NP, NC );
  auto dPhaseDens_dC     = invertLayout( phaseDens.dComp, NP, NC );
  auto dPhaseVisc_dC     = invertLayout( phaseVisc.dComp, NP, NC );
  auto dTotalDens_dC     = invertLayout( totalDens.dComp, NC);
  auto dPhaseCompFrac_dC = invertLayout( phaseCompFrac.dComp, NP, NC, NC );

  array1d<real64> compNew( NC );
  for (localIndex jc = 0; jc < NC; ++jc)
  {
    real64 const dC = perturbParameter * (composition[jc] + perturbParameter);
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compNew[ic] = composition[ic];
    }
    compNew[jc] += dC;

    // renormalize
    real64 sum = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
      sum += compNew[ic];
    for (localIndex ic = 0; ic < NC; ++ic)
      compNew[ic] /= sum;

    fluidCopy->StateUpdatePointMultiFluid( P, T, compNew, 0, 0 );
    string var = "compFrac[" + components[jc] + "]";

    checkDerivative( phaseFracCopy, phaseFrac.value, dPhaseFrac_dC[jc], dC, relTol, "phaseFrac", var, phases );
    checkDerivative( phaseDensCopy, phaseDens.value, dPhaseDens_dC[jc], dC, relTol, "phaseDens", var, phases );
    checkDerivative( phaseViscCopy, phaseVisc.value, dPhaseVisc_dC[jc], dC, relTol, "phaseVisc", var, phases );
    checkDerivative( totalDensCopy, totalDens.value, dTotalDens_dC[jc], dC, relTol, "totalDens", var );
    checkDerivative( phaseCompFracCopy, phaseCompFrac.value, dPhaseCompFrac_dC[jc], dC, relTol,
                     "phaseCompFrac", var, phases, components );
  }
}

MultiFluidBase * makeCompositionalFluid( string const & name, ManagedGroup * parent )
{
  auto fluid = parent->RegisterGroup<CompositionalMultiphaseFluid>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  compNames.resize( 4 );
  compNames[0] = "N2"; compNames[1] = "C10"; compNames[2] = "C20"; compNames[3] = "H20";

  auto & molarWgt = fluid->getReference<array1d<real64>>( MultiFluidBase::viewKeyStruct::componentMolarWeightString );
  molarWgt.resize( 4 );
  molarWgt[0] = 28e-3; molarWgt[1] = 134e-3; molarWgt[2] = 275e-3; molarWgt[3] = 18e-3;

  auto & phaseNames = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & eqnOfState = fluid->getReference<string_array>( CompositionalMultiphaseFluid::viewKeyStruct::equationsOfStateString );
  eqnOfState.resize( 2 );
  eqnOfState[0] = "PR"; eqnOfState[1] = "PR";

  auto & critPres = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalPressureString );
  critPres.resize( 4 );
  critPres[0] = 34e5; critPres[1] = 25.3e5; critPres[2] = 14.6e5; critPres[3] = 220.5e5;

  auto & critTemp = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalTemperatureString );
  critTemp.resize( 4 );
  critTemp[0] = 126.2; critTemp[1] = 622.0; critTemp[2] = 782.0; critTemp[3] = 647.0;

  auto & acFactor = fluid->getReference<array1d<real64>>( CompositionalMultiphaseFluid::viewKeyStruct::componentAcentricFactorString );
  acFactor.resize( 4 );
  acFactor[0] = 0.04; acFactor[1] = 0.443; acFactor[2] = 0.816; acFactor[3] = 0.344;

  fluid->ProcessInputFile_PostProcess();
  return fluid;
}

TEST(testMultiFluid, numericalDerivatives_compositionalFluid_molar)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  MultiFluidBase * fluid = makeCompositionalFluid( "fluid", parent.get() );
  fluid->setMassFlag( false );

  parent->Initialize( parent.get() );
  parent->FinalInitializationRecursive( parent.get() );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d<real64> comp(4);
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  testNumericalDerivatives( fluid, P, T, comp, eps, tol );
}

TEST(testMultiFluid, numericalDerivatives_compositionalFluid_mass)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  MultiFluidBase * fluid = makeCompositionalFluid( "fluid", parent.get() );
  fluid->setMassFlag( true );

  parent->Initialize( parent.get() );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d<real64> comp(4);
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-2;

  testNumericalDerivatives( fluid, P, T, comp, eps, tol );
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

#ifdef __clang__
#pragma clang diagnostic pop
#endif
