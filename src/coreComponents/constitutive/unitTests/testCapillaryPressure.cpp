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
#include "constitutive/CapillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/CapillaryPressure/VanGenuchtenCapillaryPressure.hpp"

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

void testNumericalDerivatives( CapillaryPressureBase * capPressure,
                               arraySlice1d<real64> const & saturation,
                               real64 perturbParameter,
                               real64 relTol )
{
  localIndex const NP = capPressure->numFluidPhases();

  auto const & phases = capPressure->getReference<string_array>( CapillaryPressureBase::viewKeyStruct::phaseNamesString );

  // create a clone of the capillary pressure to run updates on
  auto capPressureCopyPtr = capPressure->DeliverClone( "fluidCopy", nullptr );
  auto capPressureCopy = capPressureCopyPtr->group_cast<CapillaryPressureBase *>();

  capPressure->AllocateConstitutiveData( capPressure->getParent(), 1 );
  capPressureCopy->AllocateConstitutiveData( capPressure->getParent(), 1 );

  arraySlice1d<real64> phaseCapPressure = capPressure->getReference<array3d<real64>>( CapillaryPressureBase::viewKeyStruct::phaseCapPressureString )[0][0];
  arraySlice2d<real64> dPhaseCapPressure_dSat = capPressure->getReference<array4d<real64>>( CapillaryPressureBase::viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString )[0][0];

  arraySlice1d<real64> phaseCapPressureCopy = capPressureCopy->getReference<array3d<real64>>( CapillaryPressureBase::viewKeyStruct::phaseCapPressureString )[0][0];

  // set the fluid state to current
  capPressure->PointUpdate(saturation, 0, 0);
  
  // update saturation and check derivatives
  auto dPhaseCapPressure_dS = invertLayout( dPhaseCapPressure_dSat, NP, NP );

  array1d<real64> satNew( NP );
  for (localIndex jp = 0; jp < NP; ++jp)
  {
    real64 const dS = perturbParameter * (saturation[jp] + perturbParameter);
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      satNew[ip] = saturation[ip];
    }
    satNew[jp] += dS;

    capPressureCopy->PointUpdate( satNew, 0, 0 );
    string var = "phaseVolFrac[" + phases[jp] + "]";
    
    checkDerivative( phaseCapPressureCopy, phaseCapPressure, dPhaseCapPressure_dS[jp], dS, relTol, "phaseCapPressure", var, phases );
  }
}


CapillaryPressureBase * makeBrooksCoreyCapPressureTwoPhase( string const & name, ManagedGroup * parent )
{
  auto capPressure = parent->RegisterGroup<BrooksCoreyCapillaryPressure>( name );

  auto & phaseNames = capPressure->getReference<string_array>( CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  auto & phaseMinSat = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.05;

  auto & phaseCapPressureExpInv = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = 4;

  auto & phaseEntryPressure = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = 1;

  auto & capPressureEpsilon = capPressure->getReference<real64>( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4;
  
  capPressure->ProcessInputFile_PostProcess();
  return capPressure;
}


CapillaryPressureBase * makeBrooksCoreyCapPressureThreePhase( string const & name, ManagedGroup * parent )
{
  auto capPressure = parent->RegisterGroup<BrooksCoreyCapillaryPressure>( name );

  auto & phaseNames = capPressure->getReference<string_array>( CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "water"; phaseNames[1] = "oil"; phaseNames[2] = "gas";

  auto & phaseMinSat = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.02; phaseMinSat[2] = 0.1;

  auto & phaseCapPressureExpInv = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = -3; phaseCapPressureExpInv[2] = 2.5;

  auto & phaseEntryPressure = capPressure->getReference<array1d<real64>>( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString );
  phaseEntryPressure.resize( 3 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = -1; phaseEntryPressure[2] = 2;

  auto & capPressureEpsilon = capPressure->getReference<real64>( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-7;
  
  capPressure->ProcessInputFile_PostProcess();
  return capPressure;
}


CapillaryPressureBase * makeVanGenuchtenCapPressureTwoPhase( string const & name, ManagedGroup * parent )
{
  auto capPressure = parent->RegisterGroup<VanGenuchtenCapillaryPressure>( name );

  auto & phaseNames = capPressure->getReference<string_array>( CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & phaseMinSat = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1;

  auto & phaseCapPressureExpInv = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 0.4; phaseCapPressureExpInv[1] = 0.5;

  auto & phaseCapPressureMultiplier = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString );
  phaseCapPressureMultiplier.resize( 2 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1;

  auto & capPressureEpsilon = capPressure->getReference<real64>( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4; 

  capPressure->ProcessInputFile_PostProcess();
  return capPressure;
}

CapillaryPressureBase * makeVanGenuchtenCapPressureThreePhase( string const & name, ManagedGroup * parent )
{
  auto capPressure = parent->RegisterGroup<VanGenuchtenCapillaryPressure>( name );

  auto & phaseNames = capPressure->getReference<string_array>( CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & phaseMinSat = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1; phaseMinSat[2] = 0.;

  auto & phaseCapPressureExpInv = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 0.33; phaseCapPressureExpInv[1] = 0.4; phaseCapPressureExpInv = 0.5;

  auto & phaseCapPressureMultiplier = capPressure->getReference<array1d<real64>>( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString );
  phaseCapPressureMultiplier.resize( 3 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1; phaseCapPressureMultiplier[2] = 0.2;

  auto & capPressureEpsilon = capPressure->getReference<real64>( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4; 

  capPressure->ProcessInputFile_PostProcess();
  return capPressure;
}


TEST(testCapPressure, numericalDerivatives_brooksCoreyCapPressureTwoPhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid = makeBrooksCoreyCapPressureTwoPhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->FinalInitializationRecursive( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  real64 const start_sat = 0.2; 
  real64 const end_sat   = 0.8;
  real64 const dS = 1e-4;
  array1d<real64> sat(2);
  sat[0] = start_sat; sat[1] = 1.0-sat[0];
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}


TEST(testCapPressure, numericalDerivatives_brooksCoreyCapPressureThreePhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid = makeBrooksCoreyCapPressureThreePhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->FinalInitializationRecursive( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  real64 const start_sat = 0.2;
  real64 const end_sat   = 0.8;
  real64 const dS = 1e-4;
  array1d<real64> sat(3);
  sat[0] = start_sat; sat[1] = 0.5*(1-sat[0]); sat[2] = 1.0-sat[0]-sat[1]; 
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1-sat[0] );
    sat[2] = 1.0 - sat[0] - sat[1];
  }
}


TEST(testCapPressure, numericalDerivatives_vanGenuchtenCapPressureTwoPhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid = makeVanGenuchtenCapPressureTwoPhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->FinalInitializationRecursive( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  real64 const start_sat = 0.2;
  real64 const end_sat   = 0.8;
  real64 const dS        = 1e-4;
  array1d<real64> sat(2);
  sat[0] = start_sat; sat[1] = 1-sat[1];
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
     
}


TEST(testCapPressure, numericalDerivatives_vanGenuchtenCapPressureThreePhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid = makeVanGenuchtenCapPressureThreePhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->FinalInitializationRecursive( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  real64 const start_sat = 0.2;
  real64 const end_sat   = 0.8;
  real64 const dS        = 1e-4;
  array1d<real64> sat(2);
  sat[0] = start_sat; sat[1] = 0.5*(1-sat[0]); sat[2] = 1.0-sat[0]-sat[1]; 
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
     
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
