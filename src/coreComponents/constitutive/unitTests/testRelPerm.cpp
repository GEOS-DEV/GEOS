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
#include "constitutive/RelPerm/BrooksCoreyRelativePermeability.hpp"
#include "constitutive/RelPerm/BrooksCoreyBakerRelativePermeability.hpp"
#include "constitutive/RelPerm/VanGenuchtenBakerRelativePermeability.hpp"

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

void testNumericalDerivatives( RelativePermeabilityBase * relPerm,
                               arraySlice1d<real64> const & saturation,
                               real64 perturbParameter,
                               real64 relTol )
{
  localIndex const NP = relPerm->numFluidPhases();

  auto const & phases = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );

  // create a clone of the rel perm to run updates on
  auto relPermCopyPtr = relPerm->DeliverClone( "fluidCopy", nullptr );
  auto relPermCopy = relPermCopyPtr->group_cast<RelativePermeabilityBase *>();

  relPerm->AllocateConstitutiveData( relPerm->getParent(), 1 );
  relPermCopy->AllocateConstitutiveData( relPerm->getParent(), 1 );

  arraySlice1d<real64> phaseRelPerm = relPerm->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString )[0][0];
  arraySlice2d<real64> dPhaseRelPerm_dSat = relPerm->getReference<array4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString )[0][0];

  arraySlice1d<real64> phaseRelPermCopy = relPermCopy->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString )[0][0];

  // set the fluid state to current
  relPerm->PointUpdate( saturation, 0, 0 );

  // update saturation and check derivatives
  auto dPhaseRelPerm_dS = invertLayout( dPhaseRelPerm_dSat, NP, NP );

  array1d<real64> satNew( NP );
  for (localIndex jp = 0; jp < NP; ++jp)
  {
    real64 const dS = perturbParameter * (saturation[jp] + perturbParameter);
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      satNew[ip] = saturation[ip];
    }
    satNew[jp] += dS;

    relPermCopy->PointUpdate( satNew, 0, 0 );
    string var = "phaseVolFrac[" + phases[jp] + "]";

    checkDerivative( phaseRelPermCopy, phaseRelPerm, dPhaseRelPerm_dS[jp], dS, relTol, "phaseRelPerm", var, phases );
  }
}

RelativePermeabilityBase * makeBrooksCoreyRelPerm( string const & name, ManagedGroup * parent )
{
  auto relPerm = parent->RegisterGroup<BrooksCoreyRelativePermeability>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & phaseNames = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & phaseMinSat = relPerm->getReference<array1d<real64>>( BrooksCoreyRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.15;

  auto & phaseRelPermExp = relPerm->getReference<array1d<real64>>( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermExponentString );
  phaseRelPermExp.resize( 2 );
  phaseRelPermExp[0] = 2.0; phaseRelPermExp[1] = 2.0;

  auto & phaseRelPermMaxVal = relPerm->getReference<array1d<real64>>( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermMaxValueString );
  phaseRelPermMaxVal.resize( 2 );
  phaseRelPermMaxVal[0] = 0.8; phaseRelPermMaxVal[1] = 0.9;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeBrooksCoreyBakerRelPermTwoPhase( string const & name, ManagedGroup * parent )
{
  auto relPerm = parent->RegisterGroup<BrooksCoreyBakerRelativePermeability>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & phaseNames = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil"; 

  auto & phaseMinSat = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; 

  auto & waterOilRelPermExp = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 1.9; waterOilRelPermExp[1] = 3.95;

  auto & waterOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.8; waterOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}


RelativePermeabilityBase * makeBrooksCoreyBakerRelPermThreePhase( string const & name, ManagedGroup * parent )
{
  auto relPerm = parent->RegisterGroup<BrooksCoreyBakerRelativePermeability>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & phaseNames = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & phaseMinSat = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  auto & waterOilRelPermExp = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 2.4; waterOilRelPermExp[1] = 2.5;

  auto & waterOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.95;

  auto & gasOilRelPermExp = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentString );
  gasOilRelPermExp.resize( 2 );
  gasOilRelPermExp[0] = 1.9; gasOilRelPermExp[1] = 3.95;

  auto & gasOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase * makeVanGenuchtenBakerRelPermTwoPhase( string const & name, ManagedGroup * parent )
{
  auto relPerm = parent->RegisterGroup<VanGenuchtenBakerRelativePermeability>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & phaseNames = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; 

  auto & phaseMinSat = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.02; phaseMinSat[1] = 0.05; 

  auto & gasOilRelPermExpInv = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.7; gasOilRelPermExpInv[1] = 2.15;

  auto & gasOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.5; gasOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}


RelativePermeabilityBase * makeVanGenuchtenBakerRelPermThreePhase( string const & name, ManagedGroup * parent )
{
  auto relPerm = parent->RegisterGroup<VanGenuchtenBakerRelativePermeability>( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & phaseNames = relPerm->getReference<string_array>( RelativePermeabilityBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & phaseMinSat = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  auto & waterOilRelPermExpInv = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentInvString );
  waterOilRelPermExpInv.resize( 2 );
  waterOilRelPermExpInv[0] = 2.4; waterOilRelPermExpInv[1] = 2.5;

  auto & waterOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.95;

  auto & gasOilRelPermExpInv = relPerm->getReference<array1d<real64>>( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.9; gasOilRelPermExpInv[1] = 3.95;

  auto & gasOilRelPermMaxVal = relPerm->getReference<array1d<real64>>( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.75;

  relPerm->PostProcessInputRecursive();
  return relPerm;
}



TEST(testRelPerm, numericalDerivatives_brooksCoreyRelPerm)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * fluid = makeBrooksCoreyRelPerm( "relPerm", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  // TODO test over a range of values
  array1d<real64> sat(4);
  sat[0] = 0.7; sat[1] = 0.3;

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  testNumericalDerivatives( fluid, sat, eps, tol );
}

TEST(testRelPerm, numericalDerivatives_BrooksCoreyBakerRelPermTwoPhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * fluid = makeBrooksCoreyBakerRelPermTwoPhase( "relPerm", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-2;
  real64 const alpha = 0.4;
  array1d<real64> sat(2);
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST(testRelPerm, numericalDerivatives_BrooksCoreyBakerRelPermThreePhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * fluid = makeBrooksCoreyBakerRelPermThreePhase( "relPerm", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-2;
  real64 const alpha = 0.4;
  array1d<real64> sat(3);
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.3 *(1-sat[0]);
    sat[2] = 0.7 *(1-sat[0]);
  }
}


TEST(testRelPerm, numericalDerivatives_VanGenuchtenBakerRelPermTwoPhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * fluid = makeVanGenuchtenBakerRelPermTwoPhase( "relPerm", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-2;
  real64 const alpha = 0.4;
  array1d<real64> sat(2);
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST(testRelPerm, numericalDerivatives_VanGenuchtenBakerRelPermThreePhase)
{
  auto parent = std::make_unique<ManagedGroup>( "parent", nullptr );
  parent->resize( 1 );

  RelativePermeabilityBase * fluid = makeVanGenuchtenBakerRelPermThreePhase( "relPerm", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-2;
  real64 const alpha = 0.4;
  array1d<real64> sat(3);
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while (sat[0] <= end_sat)
  {
    testNumericalDerivatives( fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.6 *(1-sat[0]);
    sat[2] = 0.4 *(1-sat[0]);
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
