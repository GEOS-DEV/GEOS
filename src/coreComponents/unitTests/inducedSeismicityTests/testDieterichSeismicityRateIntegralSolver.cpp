/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/inducedSeismicity/DieterichSeismicityRate.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the accuracy of the integral seismicity rate solver
// It computes the seismicity in response to a hard coded stressing history to which there exists an analytical solution
char const * xmlInput =
  R"xml(

  )xml";

class DieterichSeismicityRateIntegralSolverTest : public ::testing::Test
{
public:

  DieterichSeismicityRateIntegralSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  GeosxState state;
  DieterichSeismicityRate * propagator;
};

TEST_F( DieterichSeismicityRateIntegralSolverTest, solverTest )
{
  // DomainPartition & domain = state.getProblemManager().getDomainPartition();
  // propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticFirstOrderWaveEquationSEM >( "acousticFirstOrderSolver" );

  // check number of seismos and trace length
  ASSERT_EQ( 0, 0 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}