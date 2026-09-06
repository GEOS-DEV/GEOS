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

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshManager.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/consistency/ConsistencyAdaptation.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// mixed MFD problem on an InternalMesh: the layers are then run directly with chosen parameters
static string makeInput( string const & elementType,
                         string const & xCoords,
                         string const & nx,
                         string const & cellBlocks,
                         string const & prescriptions )
{
  return GEOS_FMT( R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <InternalMesh name="mesh" elementTypes="{{ {} }}" xCoords="{{ {} }}" yCoords="{{ 0, 1 }}" zCoords="{{ 0, 1 }}"
                  nx="{{ {} }}" ny="{{ 2 }}" nz="{{ 2 }}" cellBlockNames="{{ {} }}"/>
  </Mesh>
  <Geometry>
    <Box name="thin" xMin="{{ -1e-6, -1e-6, -1e-6 }}" xMax="{{ 1.5e-3, 1.1, 1.1 }}"/>
    <Box name="east" xMin="{{ 0.9, -1e-6, -1e-6 }}" xMax="{{ 1.1, 1.1, 1.1 }}"/>
  </Geometry>
  <ElementRegions>
    <CellElementRegion name="Domain" cellBlocks="{{ * }}" materialList="{{ rock, fluid }}"/>
  </ElementRegions>
  <Solvers>
    <SinglePhaseMixedMFD name="flow" discretization="mixedMFD" targetRegions="{{ Domain }}"/>
  </Solvers>
  <NumericalMethods>
    <MixedMimetic>
      <MixedMimeticDiscretization name="mixedMFD" innerProductType="RT"/>
    </MixedMimetic>
  </NumericalMethods>
  <Constitutive>
    <CompressibleSinglePhaseFluid name="fluid" defaultDensity="1000.0" referenceDensity="1000.0"
                                  defaultViscosity="0.001" referenceViscosity="0.001" referencePressure="0.0"
                                  compressibility="0.0" viscosibility="0.0" densityModelType="exponential"/>
    <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid"
                                           porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>
    <NullModel name="nullSolid"/>
    <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.1" referencePressure="0.0" compressibility="0.0"/>
    <ConstantPermeability name="rockPerm" permeabilityComponents="{{ 1.0e-13, 1.0e-13, 1.0e-13 }}"/>
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification name="initialPressure" initialCondition="1" setNames="{{ all }}"
                        objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0e7"/>
    {}
  </FieldSpecifications>
  <Events maxTime="1.0">
    <PeriodicEvent name="solverApplications" target="/Solvers/flow"/>
  </Events>
</Problem>)xml", elementType, xCoords, nx, cellBlocks,
                   prescriptions );
}

static string const prescriptions =
  R"xml(<FieldSpecification name="thinConsistent" initialCondition="1" setNames="{ thin }"
                        objectPath="ElementRegions/Domain" fieldName="prescribedMfdFlag" scale="1.0"/>
    <FieldSpecification name="eastDiagonal" initialCondition="1" setNames="{ east }"
                        objectPath="ElementRegions/Domain" fieldName="prescribedMfdFlag" scale="0.0"/>)xml";

static string const noPrescription;

// Cartesian hexahedra: K-orthogonal, the two-point product is consistent everywhere
static string const cartesianHexa = makeInput( "C3D8", "0, 1", "10", "cb1", noPrescription );
// tetrahedra: the two-point product is not consistent in any cell
static string const tetrahedra = makeInput( "C3D4", "0, 1", "4", "cb1", noPrescription );
// a 1e-3 thick layer of 4 cells (0.2475 % of their node star) in front of 40 unit-height cells
static string const thinLayer = makeInput( "C3D8", "0, 1e-3, 1", "1, 10", "cb1, cb2", noPrescription );
static string const thinLayerPrescribed = makeInput( "C3D8", "0, 1e-3, 1", "1, 10", "cb1, cb2", prescriptions );

static void setupProblemFromXML( ProblemManager & problemManager, string const & xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult const xmlResult = xmlDocument.loadString( xmlInput );
  ASSERT_TRUE( xmlResult ) << xmlResult.description();

  Group & commandLine = problemManager.getGroup< Group >( problemManager.groupKeys.commandLine );
  commandLine.registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( MpiWrapper::commSize( MPI_COMM_GEOS ) );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( keys::ProblemManager );
  problemManager.processInputFileRecursive( xmlDocument, xmlProblemNode );

  DomainPartition & domain = problemManager.getDomainPartition();
  constitutive::ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str() );
  constitutiveManager.processInputFileRecursive( xmlDocument, topLevelNode );

  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( xmlDocument, topLevelNode );

  problemManager.problemSetup();
  problemManager.applyInitialConditions();
}

class ConsistencyAdaptationTest
{
public:

  explicit ConsistencyAdaptationTest( string const & xml )
    : m_state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {
    setupProblemFromXML( m_state.getProblemManager(), xml );
    m_regionNames.emplace_back( "Domain" );
    m_regionFilter.insert( elemManager().getRegions().getIndex( "Domain" ) );
  }

  DomainPartition & domain() { return m_state.getProblemManager().getDomainPartition(); }
  MeshLevel & mesh() { return domain().getMeshBody( 0 ).getBaseDiscretization(); }
  ElementRegionManager & elemManager() { return mesh().getElemManager(); }

  /// run the layers and reduce the counts over the ranks
  ConsistencyAdaptation::Report classify( ConsistencyAdaptation::Parameters const & params )
  {
    ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > const permeability =
      elemManager().constructMaterialArrayViewAccessor< constitutive::PermeabilityBase, real64, 3 >( fields::permeability::permeability::key() );
    ConsistencyAdaptation::Report r = ConsistencyAdaptation::classify( mesh(), m_regionNames, m_regionFilter.toViewConst(),
                                                                       permeability.toNestedViewConst(), params, domain().getNeighbors() );
    r.numCells = MpiWrapper::sum( r.numCells );
    r.numConsistent = MpiWrapper::sum( r.numConsistent );
    r.numPrescribed0 = MpiWrapper::sum( r.numPrescribed0 );
    r.numPrescribed1 = MpiWrapper::sum( r.numPrescribed1 );
    r.numDegenerate = MpiWrapper::sum( r.numDegenerate );
    r.numRejected = MpiWrapper::sum( r.numRejected );
    return r;
  }

  /// number of owned cells with eta = value
  localIndex countFlag( integer const value )
  {
    localIndex count = 0;
    elemManager().forElementSubRegions< ElementSubRegionBase >( m_regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const mfdFlag = subRegion.getField< fields::mixedMimetic::mfdFlag >();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        count += ( ghostRank[ei] < 0 && mfdFlag[ei] == value ) ? 1 : 0;
      }
    } );
    return MpiWrapper::sum( count );
  }

  /// number of owned cells whose degeneracy indicator is below the percentage
  localIndex countDegeneracyBelow( real64 const percent )
  {
    localIndex count = 0;
    elemManager().forElementSubRegions< ElementSubRegionBase >( m_regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const indicator = subRegion.getField< fields::mixedMimetic::degeneracyIndicator >();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        count += ( ghostRank[ei] < 0 && indicator[ei] < percent ) ? 1 : 0;
      }
    } );
    return MpiWrapper::sum( count );
  }

  /// check label_f = max of eta over the cells of f, and return the number of owned faces with label 1
  localIndex checkFaceLabels()
  {
    FaceManager const & faceManager = mesh().getFaceManager();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    arrayView1d< integer const > const label = faceManager.getField< fields::mixedMimetic::faceStencilLabel >();
    arrayView1d< integer const > const ghostRank = faceManager.ghostRank();
    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const mfdFlag =
      elemManager().constructArrayViewAccessor< integer, 1 >( fields::mixedMimetic::mfdFlag::key() );

    localIndex count = 0;
    for( localIndex kf = 0; kf < faceManager.size(); ++kf )
    {
      integer expected = 0;
      for( localIndex k = 0; k < elemRegionList.size( 1 ); ++k )
      {
        localIndex const er = elemRegionList( kf, k );
        localIndex const esr = elemSubRegionList( kf, k );
        localIndex const ei = elemList( kf, k );
        if( er >= 0 && esr >= 0 && ei >= 0 && m_regionFilter.contains( er ) )
        {
          expected = std::max( expected, mfdFlag[er][esr][ei] );
        }
      }
      EXPECT_EQ( label[kf], expected ) << "face " << kf;
      count += ( ghostRank[kf] < 0 && label[kf] == 1 ) ? 1 : 0;
    }
    return MpiWrapper::sum( count );
  }

private:

  GeosxState m_state;
  string_array m_regionNames;
  SortedArray< localIndex > m_regionFilter;
};

TEST( ConsistencyAdaptation, ConsistencyLayer_CartesianHexahedra )
{
  ConsistencyAdaptationTest test( cartesianHexa );
  ConsistencyAdaptation::Parameters params;

  // the two-point product is consistent on K-orthogonal cells: eta = 0 everywhere (up to round-off)
  params.consistencyTolerance = 1e-10;
  ConsistencyAdaptation::Report r = test.classify( params );
  EXPECT_EQ( r.numCells, 40 );
  EXPECT_EQ( r.numConsistent, 0 );
  EXPECT_EQ( test.countFlag( 1 ), 0 );
  EXPECT_EQ( test.checkFaceLabels(), 0 );

  // without the layer the selected product is used everywhere
  params.adaptiveConsistency = false;
  r = test.classify( params );
  EXPECT_EQ( r.numConsistent, 40 );
  EXPECT_EQ( test.countFlag( 1 ), 40 );
  EXPECT_EQ( test.checkFaceLabels(), 164 );

  // a TPFA discretization condenses every face whatever eta
  params.effectiveTpfa = true;
  test.classify( params );
  EXPECT_EQ( test.countFlag( 1 ), 40 );
  FaceManager const & faceManager = test.mesh().getFaceManager();
  arrayView1d< integer const > const label = faceManager.getField< fields::mixedMimetic::faceStencilLabel >();
  for( localIndex kf = 0; kf < faceManager.size(); ++kf )
  {
    EXPECT_EQ( label[kf], 0 );
  }
}

TEST( ConsistencyAdaptation, ConsistencyLayer_Tetrahedra )
{
  ConsistencyAdaptationTest test( tetrahedra );
  ConsistencyAdaptation::Parameters params;

  // the two-point product is not consistent on any simplex
  params.consistencyTolerance = 1e-12;
  ConsistencyAdaptation::Report r = test.classify( params );
  EXPECT_EQ( r.numCells, 96 );
  EXPECT_EQ( r.numConsistent, 96 );
  EXPECT_EQ( test.countFlag( 1 ), 96 );
  EXPECT_EQ( test.checkFaceLabels(), MpiWrapper::sum( [&]
  {
    arrayView1d< integer const > const ghostRank = test.mesh().getFaceManager().ghostRank();
    localIndex owned = 0;
    for( localIndex kf = 0; kf < ghostRank.size(); ++kf )
    {
      owned += ( ghostRank[kf] < 0 ) ? 1 : 0;
    }
    return owned;
  }() ) );

  // the number of consistent cells is non-increasing with the tolerance
  localIndex previous = 96;
  for( real64 const tolerance : { 1e-3, 1e-1, 1.0, 1e+20 } )
  {
    params.consistencyTolerance = tolerance;
    r = test.classify( params );
    EXPECT_LE( r.numConsistent, previous ) << "tolerance = " << tolerance;
    EXPECT_EQ( r.numConsistent, test.countFlag( 1 ) );
    test.checkFaceLabels();
    previous = r.numConsistent;
  }
  EXPECT_EQ( previous, 0 );
}

TEST( ConsistencyAdaptation, DegeneracyLayer_ThinLayer )
{
  ConsistencyAdaptationTest test( thinLayer );
  ConsistencyAdaptation::Parameters params;
  params.adaptiveConsistency = false;

  // the 4 thin cells are 0.2475 % of their node star
  params.degeneracyTolerance = 0.0;
  ConsistencyAdaptation::Report r = test.classify( params );
  EXPECT_EQ( r.numCells, 44 );
  EXPECT_EQ( r.numDegenerate, 0 );
  EXPECT_EQ( test.countFlag( 1 ), 44 );
  EXPECT_EQ( test.countDegeneracyBelow( 0.2 ), 0 );
  EXPECT_EQ( test.countDegeneracyBelow( 0.3 ), 4 );
  EXPECT_EQ( test.countDegeneracyBelow( 100.0 ), 44 );

  params.degeneracyTolerance = 0.2;
  r = test.classify( params );
  EXPECT_EQ( r.numDegenerate, 0 );
  EXPECT_EQ( test.countFlag( 1 ), 44 );

  // the 16 faces of the thin block not shared with a thick cell are condensed
  params.degeneracyTolerance = 1.0;
  r = test.classify( params );
  EXPECT_EQ( r.numDegenerate, 4 );
  EXPECT_EQ( r.numRejected, 0 );
  EXPECT_EQ( test.countFlag( 1 ), 40 );
  EXPECT_EQ( test.checkFaceLabels(), 180 - 16 );
}

TEST( ConsistencyAdaptation, Prescription_ThinLayer )
{
  ConsistencyAdaptationTest test( thinLayerPrescribed );
  ConsistencyAdaptation::Parameters params;

  // the prescription overrides the consistency layer (eta = 0 everywhere on this mesh)
  params.consistencyTolerance = 1e-10;
  params.degeneracyTolerance = 0.0;
  ConsistencyAdaptation::Report r = test.classify( params );
  EXPECT_EQ( r.numConsistent, 0 );
  EXPECT_EQ( r.numPrescribed1, 4 );
  EXPECT_EQ( r.numPrescribed0, 4 );
  EXPECT_EQ( r.numDegenerate, 0 );
  EXPECT_EQ( test.countFlag( 1 ), 4 );
  EXPECT_EQ( test.checkFaceLabels(), 20 );

  // and the selected product everywhere: the 4 east cells are prescribed the diagonal product,
  // the 16 faces of the east block not shared with the previous column are condensed
  params.adaptiveConsistency = false;
  r = test.classify( params );
  EXPECT_EQ( test.countFlag( 0 ), 4 );
  EXPECT_EQ( test.checkFaceLabels(), 180 - 16 );

  // only the degeneracy layer rejects a prescribed consistent product
  params.adaptiveConsistency = true;
  params.degeneracyTolerance = 1.0;
  r = test.classify( params );
  EXPECT_EQ( r.numPrescribed1, 4 );
  EXPECT_EQ( r.numDegenerate, 4 );
  EXPECT_EQ( r.numRejected, 4 );
  EXPECT_EQ( test.countFlag( 1 ), 0 );
  EXPECT_EQ( test.checkFaceLabels(), 0 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
