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

bool testNumericalDerivatives( MultiFluidBase * fluid,
                               real64 P,
                               real64 T,
                               array1d<real64> const & composition,
                               real64 perturbParameter = 1e-6,
                               real64 relTol = 1e-4 )
{
  localIndex const NC = fluid->numFluidComponents();
  localIndex const NP = fluid->numFluidPhases();

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  // create a clone of the fluid to run updates on
  auto fluidCopy = fluid->DeliverClone( "fluidCopy", nullptr );

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

  // extract data views (without derivatives) from the copy

  // set the fluid state to current
  fluid->StateUpdatePointMultiphaseFluid( P, T, composition.data(), 0, 0 );

  bool result = true;

  // update pressure and check derivatives
  {
    real64 const dP = perturbParameter * std::fmax( P, 1.0 );
    fluidCopy->StateUpdatePointMultiphaseFluid( P + dP, T, composition.data(), 0, 0 );

    // check phase fraction derivatives
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const numDeriv = (phaseFracCopy[ip] - phaseFrac.value[ip]) / dP;
      real64 const anlDeriv = phaseFrac.dPres[ip];
      real64 const delta = std::fabs( anlDeriv - numDeriv );
      real64 const value = std::fmax( std::fabs(anlDeriv), std::fabs(numDeriv) );
      if (value > 0 && delta / value > relTol)
      {
        GEOS_LOG( "d(phaseFrac[" << phases[ip] << "])/dP mismatch in rel tol: " << delta/value << " > " << relTol );
        result = false;
      }
    }

    // check phase density derivatives
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const numDeriv = (phaseDensCopy[ip] - phaseDens.value[ip]) / dP;
      real64 const anlDeriv = phaseDens.dPres[ip];
      real64 const delta = std::fabs( anlDeriv - numDeriv );
      real64 const value = std::fmax( std::fabs(anlDeriv), std::fabs(numDeriv) );
      if (value > 0 && delta / value > relTol)
      {
        GEOS_LOG( "d(phaseDens[" << phases[ip] << "])/dP mismatch in rel tol: " << delta/value << " > " << relTol );
        result = false;
      }
    }

    // check phase viscosity derivatives
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const numDeriv = (phaseViscCopy[ip] - phaseVisc.value[ip]) / dP;
      real64 const anlDeriv = phaseVisc.dPres[ip];
      real64 const delta = std::fabs( anlDeriv - numDeriv );
      real64 const value = std::fmax( std::fabs(anlDeriv), std::fabs(numDeriv) );
      if (value > 0 && delta / value > relTol)
      {
        GEOS_LOG( "d(phaseVisc[" << phases[ip] << "])/dP mismatch in rel tol: " << delta/value << " > " << relTol );
        result = false;
      }
    }

    // check phase component fractions derivatives
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const numDeriv = (phaseCompFracCopy[ip][ic] - phaseCompFrac.value[ip][ic]) / dP;
        real64 const anlDeriv = phaseCompFrac.dPres[ip][ic];
        real64 const delta = std::fabs(anlDeriv - numDeriv);
        real64 const value = std::fmax(std::fabs(anlDeriv), std::fabs(numDeriv));
        if (value > 0 && delta / value > relTol)
        {
          GEOS_LOG( "d(phaseCompFrac[" << phases[ip] << "][" << components[ic] << "])/dP mismatch in rel tol: "
                    << delta / value << " > " << relTol );
          result = false;
        }
      }
    }

    // check total density derivatives
    {
      real64 const numDeriv = (totalDensCopy - totalDens.value) / dP;
      real64 const anlDeriv = totalDens.dPres;
      real64 const delta = std::fabs( anlDeriv - numDeriv );
      real64 const value = std::fmax( std::fabs(anlDeriv), std::fabs(numDeriv) );
      if (value > 0 && delta / value > relTol)
      {
        GEOS_LOG( "d(totalDens)/dP mismatch in rel tol: " << delta/value << " > " << relTol );
        result = false;
      }
    }
  };

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

TEST(testMultiFluid, numericalDerivativesCompMutliFluid)
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

  bool derivsAreOk = testNumericalDerivatives( fluid, P, T, comp, 1e-8, 1e-5 );

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
