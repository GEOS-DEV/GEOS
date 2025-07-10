
#include <gtest/gtest.h>
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/gravity/GravityFE.hpp"
#include "mesh/DomainPartition.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"


#include "common/DataTypes.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;

using namespace geos::finiteElement;

CommandLineOptions g_commandLineOptions;

// Reusable adjoint test function
template< typename ForwardOp, typename AdjointOp >
std::pair< bool, double > adjointTest(
  ForwardOp A,
  AdjointOp AT,
  std::vector< double > const & x,
  std::vector< double > const & y,
  double tol = 1e-13,
  bool verbose = true )
{
  auto dot = []( std::vector< double > const & a, std::vector< double > const & b ) -> double {
    return std::inner_product( a.begin(), a.end(), b.begin(), 0.0 );
  };

  auto norm = [&]( std::vector< double > const & v ) -> double {
    return std::sqrt( dot( v, v ));
  };

  std::vector< double > Ax = A( x );
  std::vector< double > ATy = AT( y );

  double IP1 = dot( Ax, y );
  double IP2 = dot( x, ATy );

  double error = std::abs( IP1 - IP2 ) / std::max( norm( Ax ) * norm( y ), norm( ATy ) * norm( x ));
  bool passed = error < tol;

  if( verbose )
  {
    std::cout << "\n=== Adjoint Test Summary ===\n";
    std::cout << "<Ax, y>     = " << IP1 << "\n";
    std::cout << "<x, A^T y>  = " << IP2 << "\n";
    std::cout << "Passed      = " << std::boolalpha << passed << "\n";
    std::cout << "Error       = " << error << "\n";
  }

  return { passed, error };
}

TEST( GravityFEKernelTest, AdjointConsistencyWithGravityFE )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ));

  const char * xmlInput =
    "<?xml version=\"1.0\" ?>"
    "<Problem>"
    "  <Solvers>"
    "    <GravityFE name='gravity' discretization='FE1' targetRegions='{ region }' stationCoordinates='{ { 0.5, 0.5, 2.0 }, { 1.5, 0.5, 2.0 }, { 2.5, 0.5, 2.0 } }'/>"
    "  </Solvers>"
    "  <Mesh>"
    "    <InternalMesh name='mesh' elementTypes='{ C3D8 }' xCoords='{ 10000, 14000 }' yCoords='{ 8000, 11000 }' zCoords='{ 2000, 2400 }'"
    "      nx='{ 12 }' ny='{ 12 }' nz='{ 12 }' cellBlockNames='{ cellBlock }'/>"
    "  </Mesh>"
    "  <Events maxTime='1.00001e6'>"
    "    <SoloEvent name='solverApplicationsGravity' targetTime='1e6' target='/Solvers/gravity'/>"
    "  </Events>"
    "  <NumericalMethods>"
    "    <FiniteElements>"
    "      <FiniteElementSpace name='FE1' order='1'/>"
    "    </FiniteElements>"
    "  </NumericalMethods>"
    "  <ElementRegions>"
    "    <CellElementRegion name='region' cellBlocks='{ cellBlock }' materialList='{ nullModel }'/>"
    "  </ElementRegions>"
    "  <Constitutive>"
    "    <NullModel name='nullModel'/>"
    "  </Constitutive>"
    "  <FieldSpecifications>"
    "    <FieldSpecification name='mediumDensity' initialCondition='1' setNames='{ all }' objectPath='ElementRegions' fieldName='mediumDensity' scale='900'/>"
    "  </FieldSpecifications>"
    "</Problem>";

  setupProblemFromXML( state.getProblemManager(), xmlInput );

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  GravityFE & gravityFE = state.getProblemManager().getPhysicsSolverManager().getGroup< GravityFE >( "gravity" );

  // Set up synthetic station and residue
  array2d< real64 > & stations = gravityFE.getReference< array2d< real64 > >( "stationCoordinates" );
  stations.resize( 3, 3 );
  stations( 0, 0 ) = 12557.71; stations( 0, 1 ) = 8075.03; stations( 0, 2 ) = 2527.50;
  stations( 1, 0 ) = 10892.84; stations( 1, 1 ) = 10209.41; stations( 1, 2 ) = 2567.67;
  stations( 2, 0 ) = 13568.72; stations( 2, 1 ) = 8260.82; stations( 2, 2 ) = 2542.19;



  array1d< real64 > & residueField = gravityFE.getReference< array1d< real64 > >( "residue" );

  residueField.resize( 3 );
  residueField[0] = 1.0;
  residueField[1] = 1.0;
  residueField[2] = 1.0;


  // Forward operator
  auto A = [&]( std::vector< double > const & xVec ) -> std::vector< double >
  {
    gravityFE.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( std::string const &,
                                                                           MeshLevel & mesh,
                                                                           string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & elementSubRegion )
      {
        arrayView1d< real64 > const density = elementSubRegion.getReference< array1d< real64 > >( fields::MediumDensity::key());
        for( localIndex i = 0; i < density.size(); ++i )
        {
          const_cast< real64 & >(density[i]) = xVec[i];
        }
      } );
    } );



    gravityFE.explicitStepModeling( 0.0, 1.0, 0, domain );

    arrayView1d< real64 const > gz = gravityFE.getReference< array1d< real64 > >( "gzAtStations" ).toViewConst();
    return std::vector< double >( gz.begin(), gz.end());
  };



// Adjoint operator
  auto AT = [&]( std::vector< double > const & yVec ) -> std::vector< double >
  {
    // Set the residue field from yVec
    array1d< real64 > & residueInner = gravityFE.getReference< array1d< real64 > >( "residue" );
    for( localIndex i = 0; i < residueInner.size(); ++i )
      residueInner[i] = yVec[i ];

    // Run the adjoint step
    gravityFE.explicitStepAdjoint( 0.0, 1.0, 0, domain );

    // Extract adjoint field using the same loop structure
    std::vector< double > adjointValues;
    gravityFE.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( std::string const &,
                                                                           MeshLevel & mesh,
                                                                           string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & elementSubRegion )
      {
        arrayView1d< real64 > const adjoint = elementSubRegion.getReference< array1d< real64 > >( fields::Adjoint::key());
        for( localIndex i = 0; i < adjoint.size(); ++i )
        {
          adjointValues.push_back( adjoint[i] );
        }
      } );
    } );

    return adjointValues;
  };



  std::vector< double > x( 12*12*12, 1.0 ); // synthetic density

  std::vector< double > y( 3, 2.0 );  // synthetic residue



  auto [passed, error] = adjointTest( A, AT, x, y );
  EXPECT_TRUE( passed );
  EXPECT_LT( error, 1e-13 );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
