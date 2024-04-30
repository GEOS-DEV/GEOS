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
#include "physicsSolvers/inducedSeismicity/SeismicityRate.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

static constexpr real64 relTol = 1.0e-2;

// This unit test checks the accuracy of the integral seismicity rate solver
// It computes the seismicity in response to a hard coded stressing history to which there exists an analytical solution
char const * xmlInput =
  R"xml(
  <?xml version="1.0" ?>
  <Problem>
    <Solvers>
      <SeismicityRate
        name="dieterichSR"
        targetRegions="{ Domain }"
        directEffect="0.01"
        backgroundStressingRate="3.171e-5">
      </SeismicityRate>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 1e3 }"
        yCoords="{ 0, 1e3 }"
        zCoords="{ 0, 1e3 }"
        nx="{ 1 }"
        ny="{ 1 }"
        nz="{ 1 }"
        cellBlockNames="{ cb1 }"/>
    </Mesh>
    <Events
      maxTime="360000.0">
      <PeriodicEvent
        name="solverApplications"
        forceDt="3600.0"
        target="/Solvers/dieterichSR"/>
    </Events>
    <ElementRegions>
      <CellElementRegion
        name="Domain"
        cellBlocks="{ cb1 }"
        materialList="{}"/>
    </ElementRegions>
  </Problem>
  )xml";

class SeismicityRateIntegralSolverTest : public ::testing::Test
{
public:

  SeismicityRateIntegralSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 3600.0;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  // Parameters of model and analytical solution
  static real64 constexpr cTau = 1e-3;
  static real64 constexpr aSigma = 1e6;
  static real64 constexpr t_a = 1e6/3.171e-5;
  static real64 constexpr iniSig = -100e6;
  static real64 constexpr iniTau = 60e6;

  GeosxState state;
  SeismicityRate * propagator;
};

TEST_F( SeismicityRateIntegralSolverTest, solverTest )
{
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< SeismicityRate >( "dieterichSR" );
  real64 time_n = time;

  // run for 100 hours (100 steps)
  for( int i=0; i<100; i++ )
  {

    // Loop through CellElementSubRegions of domain to pass the subRegion as input to seismciity rate solver
    domain.forMeshBodies( [&] ( MeshBody & meshBody )
    {
      meshBody.forMeshLevels( [&] ( MeshLevel & mesh )
      {
        ElementRegionManager & elemManager = mesh.getElemManager();
        for( localIndex er = 0; er < elemManager.numRegions(); ++er )
        {
          ElementRegionBase & elemRegion = elemManager.getRegion( er );
          elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const, CellElementSubRegion & subRegion )
          {

            // Initialize shear and normal stresses acting on fault
            // (field variables in seismicity rate solver) as specified in XML file
            if( i==0 )
            {
              arrayView1d< real64 > const tempSigIni = subRegion.getField< fields::inducedSeismicity::initialProjectedNormalTraction >();
              tempSigIni.setValues< parallelHostPolicy >( iniSig );
              arrayView1d< real64 > const tempSig_n = subRegion.getField< fields::inducedSeismicity::projectedNormalTraction_n >();
              tempSig_n.setValues< parallelHostPolicy >( iniSig );
              arrayView1d< real64 > const tempSig = subRegion.getField< fields::inducedSeismicity::projectedNormalTraction >();
              tempSig.setValues< parallelHostPolicy >( iniSig );

              arrayView1d< real64 > const tempTauIni = subRegion.getField< fields::inducedSeismicity::initialProjectedShearTraction >();
              tempTauIni.setValues< parallelHostPolicy >( iniTau );
              arrayView1d< real64 > const tempTau_n = subRegion.getField< fields::inducedSeismicity::projectedShearTraction_n >();
              tempTau_n.setValues< parallelHostPolicy >( iniTau );
              arrayView1d< real64 > const tempTau = subRegion.getField< fields::inducedSeismicity::projectedShearTraction >();
              tempTau.setValues< parallelHostPolicy >( iniTau );
            }

            // Retrieve shear stress fields from seismicity rate solver, to be hard coded
            arrayView1d< real64 const > const tau_i = subRegion.getField< geos::fields::inducedSeismicity::initialProjectedShearTraction >();
            arrayView1d< real64 > const tau = subRegion.getField< geos::fields::inducedSeismicity::projectedShearTraction >();
            arrayView1d< real64 > const tau_n = subRegion.getField< geos::fields::inducedSeismicity::projectedShearTraction_n >();

            // Shear stress history will follow a log history,
            // to which one can compute an analytical solution for the seismicity rate
            real64 curTau = aSigma*std::log( cTau*(time_n+dt)+1 ) + tau_i[0];
            real64 curTau_n = aSigma*std::log( cTau*(time_n)+1 ) + tau_i[0];
            tau.setValues< parallelHostPolicy >( curTau );
            tau_n.setValues< parallelHostPolicy >( curTau_n );

            // Call seismicity rate solver with the modified shear stress history
            propagator->integralSolverStep( time_n, dt, subRegion );

            // retrive seismicity rate from solver
            arrayView1d< real64 const > const R = subRegion.getField< geos::fields::inducedSeismicity::seismicityRate >();

            // compute analytical solution to the hard coded shear stress history
            real64 K = std::exp( curTau/aSigma - tau_i[0]/aSigma );
            real64 denom = (cTau*(time_n+dt - t_a) + 1)*std::exp((time_n+dt)/t_a ) + cTau*t_a;
            real64 Ranly = K/denom;
            
            // Check relative error against analytical solution
            checkRelativeError( R[0], Ranly, relTol );
          } );
        }
      } );
    } );
    time_n += dt;
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
