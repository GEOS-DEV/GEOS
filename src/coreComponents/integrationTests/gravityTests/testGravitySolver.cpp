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

#include <gtest/gtest.h>

#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp" // For setupProblemFromXML
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/gravity/GravityFE.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;


// Shared constants
constexpr std::array< real64, 6 > PRISM_BOUNDS = { 10000, 14000, 8000, 11000, 2000, 2200 };
constexpr real64 DX = 500.0;
constexpr real64 DY = 500.0;
constexpr real64 Z_STATION = 10.0;
constexpr localIndex NX = 41;
constexpr localIndex NY = 31;
constexpr real64 DENSITY_CONTRAST = 900.0;


/**
 * @brief Analytic solution for gravity anomaly (gz) due to a rectangular prism.
 *
 * @param station Coordinates of the observation point (x, y, z)
 * @param prism_bounds Bounds of the prism: {x0, x1, y0, y1, z0, z1}
 * @param density_contrast Density contrast of the prism (kg/m^3)
 * @return double Gravity anomaly in m/s^2
 */
double analyticPrismGz( const std::array< double, 3 > & station,
                        const std::array< double, 6 > & prism_bounds,
                        double density_contrast )
{
  const double x = station[0];


  const double y = station[1];
  const double z = station[2];

  const double x0 = prism_bounds[0];
  const double x1 = prism_bounds[1];
  const double y0 = prism_bounds[2];
  const double y1 = prism_bounds[3];
  const double z0 = prism_bounds[4];
  const double z1 = prism_bounds[5];

  const double epsilon = 1e-16;
  double val = 0.0;

  for( int i = 0; i < 2; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      for( int k = 0; k < 2; ++k )
      {
        double xi = (i == 0) ? x0 : x1;
        double yj = (j == 0) ? y0 : y1;
        double zk = (k == 0) ? z0 : z1;

        double dx = xi - x;
        double dy = yj - y;
        double dz = zk - z;

        double r = std::sqrt( dx * dx + dy * dy + dz * dz );

        if( r > epsilon )
        {
          double sign = ((i + j + k) % 2 == 0) ? 1.0 : -1.0;
          double arg1 = std::atan2( dx * dy, dz * r + epsilon );
          double arg2 = std::log( r + dy + epsilon );
          double arg3 = std::log( r + dx + epsilon );

          val += sign * (dz * arg1 - dx * arg2 - dy * arg3);
        }
      }
    }
  }

  return -GRAVITATIONAL_CONSTANT * density_contrast * val;
}



class GravitySolverTest : public ::testing::Test
{
public:
  GravitySolverTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ))
  {}

protected:
  void SetUp() override
  {
    stations.clear();
    std::ostringstream coords;
    coords << "{ ";
    for( localIndex j = 0; j < NY; ++j )
    {
      for( localIndex i = 0; i < NX; ++i )
      {
        std::array< real64, 3 > station = { i * DX, j * DY, Z_STATION };
        stations.push_back( station );
        if( i != 0 || j != 0 ) coords << ", ";
        coords << "{ " << station[0] << ", " << station[1] << ", " << station[2] << " }";
      }
    }
    coords << " }";

    std::ostringstream xml;
    xml << "<?xml version=\"1.0\" ?>\n"
        << "<Problem>\n"
        << "  <Solvers>\n"
        << "    <GravityFE name=\"gravity\" discretization=\"FE1\" targetRegions=\"{ region }\"\n"
        << "      mode=\"modeling\" stationCoordinates=\"" << coords.str() <<"\"/>\n"
        << "  </Solvers>\n"
        << "  <Mesh>\n"
        << "    <InternalMesh name=\"mesh1\" elementTypes=\"{ C3D8 }\"\n"
        << "      xCoords=\"{ " << PRISM_BOUNDS[0] << ", " << PRISM_BOUNDS[1] << " }\"\n"
        << "      yCoords=\"{ " << PRISM_BOUNDS[2] << ", " << PRISM_BOUNDS[3] << " }\"\n"
        << "      zCoords=\"{ " << PRISM_BOUNDS[4] << ", " << PRISM_BOUNDS[5] << " }\"\n"
        << "      nx=\"{ 12 }\" ny=\"{ 12 }\" nz=\"{ 12 }\" cellBlockNames=\"{ cellBlock }\"/>\n"
        << "  </Mesh>\n"
        << "  <Events maxTime=\"1.00001e6\">\n"
        << "    <SoloEvent name=\"solverApplicationsGravity\" targetTime=\"1e6\" target=\"/Solvers/gravity\"/>\n"
        << "  </Events>\n"
        << "  <NumericalMethods>\n"
        << "    <FiniteElements>\n"
        << "      <FiniteElementSpace name=\"FE1\" order=\"1\"/>\n"
        << "    </FiniteElements>\n"
        << "  </NumericalMethods>\n"
        << "  <ElementRegions>\n"
        << "    <CellElementRegion name=\"region\" cellBlocks=\"{ cellBlock }\" materialList=\"{ nullModel }\"/>\n"
        << "  </ElementRegions>\n"
        << "  <Constitutive>\n"
        << "    <NullModel name=\"nullModel\"/>\n"
        << "  </Constitutive>\n"
        << "  <FieldSpecifications>\n"
        << "    <FieldSpecification name=\"mediumDensity\" initialCondition=\"1\" setNames=\"{ all }\"\n"
        << "      objectPath=\"ElementRegions\" fieldName=\"mediumDensity\" scale=\"" << DENSITY_CONTRAST << "\"/>\n"
        << "  </FieldSpecifications>\n"
        << "</Problem>\n";

    setupProblemFromXML( state.getProblemManager(), xml.str().c_str());
  }

  GeosxState state;
  std::vector< std::array< real64, 3 > > stations;
  GravityFE * gravitySolver;
};

TEST_F( GravitySolverTest, RectangularPrismAnomaly )
{
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  gravitySolver = &state.getProblemManager().getPhysicsSolverManager().getGroup< GravityFE >( "gravity" );

  gravitySolver->explicitStepModeling( 0.0, 1.0e6, 0, domain );

  auto const & gz_computed = gravitySolver->getReference< array1d< real64 > >( "gzAtStations" );

  EXPECT_EQ( gz_computed.size(), stations.size());

  for( localIndex i = 0; i < gz_computed.size(); ++i )
  {
    real64 gz_expected = analyticPrismGz( stations[i], PRISM_BOUNDS, DENSITY_CONTRAST );
    EXPECT_NEAR( gz_computed[i], gz_expected, 1e-6 )
      << "Station " << i << " mismatch: computed=" << gz_computed[i]
      << ", expected=" << gz_expected;
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
