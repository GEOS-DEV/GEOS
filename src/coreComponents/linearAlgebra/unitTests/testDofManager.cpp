/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testDofManager.cpp
 * @brief This test file is part of the ctest suite and tests the DofManager functionality.
 */

#include "gtest/gtest.h"

#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"

#include "testDofManagerUtils.hpp"

#include <memory>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::dataRepository;

char const * xmlInput =
"<Problem>"
"  <Mesh>"
"    <InternalMesh name=\"mesh1\""
"                  elementTypes=\"{C3D8}\""
"                  xCoords=\"{0, 1, 2, 3, 4}\""
"                  yCoords=\"{0, 1}\""
"                  zCoords=\"{0, 1}\""
"                  nx=\"{4, 4, 4, 4}\""
"                  ny=\"{4}\""
"                  nz=\"{5}\""
"                  cellBlockNames=\"{block1, block2, block3, block4}\"/>"
"  </Mesh>"
"  <ElementRegions>"
"    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{}\" />"
"    <CellElementRegion name=\"region2\" cellBlocks=\"{block2}\" materialList=\"{}\" />"
"    <CellElementRegion name=\"region3\" cellBlocks=\"{block3}\" materialList=\"{}\" />"
"    <CellElementRegion name=\"region4\" cellBlocks=\"{block4}\" materialList=\"{}\" />"
"  </ElementRegions>"
"</Problem>";

/**
 * @brief Base class for all DofManager test fixtures.
 */
class DofManagerTestBase : public ::testing::Test
{
public:

  DofManagerTestBase() :
    ::testing::Test(),
    problemManager( std::make_unique<ProblemManager>( "Problem", nullptr ) ),
    dofManager( "test" )
  {
    setupProblem( problemManager.get(), xmlInput );
    mesh = problemManager->getDomainPartition()->getMeshBody( 0 )->getMeshLevel( 0 );
    dofManager.setMesh( problemManager->getDomainPartition(), 0, 0 );
  }

protected:

  std::unique_ptr<ProblemManager> const problemManager;
  MeshLevel * mesh;
  DofManager dofManager;
};

/**
 * @brief Check that all local objects have dofs assigned. Collect all local dofs into a list.
 * @param [in] mesh the mesh
 * @param [in] dofIndexKey key used to look up dof index array
 * @param [in] regions list of target region names
 * @param [out] dofNumbers the list of dof indices
 */
template< DofManager::Location LOC >
void checkLocalDofNumbers( MeshLevel const * const mesh,
                           string const & dofIndexKey,
                           string_array const & regions,
                           array1d< globalIndex > & dofNumbers )
{
  ObjectManagerBase const * const manager = mesh->GetGroup<ObjectManagerBase>( testMeshHelper<LOC>::managerKey );
  arrayView1d< globalIndex const > dofIndex = manager->getReference< array1d<globalIndex> >( dofIndexKey );

  forLocalObjects<LOC>( mesh, regions, [&]( localIndex const idx )
  {
    SCOPED_TRACE( "idx = " + std::to_string( idx ) );
    EXPECT_GE( dofIndex[idx], 0 );
    dofNumbers.push_back( dofIndex[idx] );
  } );
}

/**
 * @brief Collect all local dofs into a list.
 * @param [in] mesh the mesh
 * @param [in] dofIndexKey key used to look up dof index array
 * @param [in] regions list of target region names
 * @param [out] dofNumbers the list of dof indices
 */
template<>
void checkLocalDofNumbers< DofManager::Location::Elem >( MeshLevel const * const mesh,
                                                         string const & dofIndexKey,
                                                         string_array const & regions,
                                                         array1d<globalIndex> & dofNumbers )
{
  // make a list of regions
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  auto const dofNumber = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofIndexKey );

  forLocalObjects<DofManager::Location::Elem>( mesh, regions, [&]( auto const idx )
  {
    globalIndex const dofIndex = dofNumber[idx[0]][idx[1]][idx[2]];
    EXPECT_GE( dofIndex, 0 );
    dofNumbers.push_back( dofIndex );
  } );
}

/**
 * @brief Checks that dof numbers are unique
 * @param [in] dofNumbers list of dof numbers, must be sorted
 */
void checkUniqueness( arrayView1d< globalIndex const > const & dofNumbers )
{
  EXPECT_TRUE( std::adjacent_find( dofNumbers.begin(), dofNumbers.end() ) == dofNumbers.end() );
}

/**
 * @brief Checks that dof numbers have a given stride between them.
 * @param [in] dofNumbers list of dof numbers, must be sorted
 * @param [in] stride expected stride between dof numbers
 */
void checkStride( arrayView1d< globalIndex const > const & dofNumbers, globalIndex const stride )
{
  // Check dof stride is equal to numComp
  for( localIndex i = 1; i < dofNumbers.size(); ++i )
  {
    SCOPED_TRACE( "i = " + std::to_string( i ) );
    EXPECT_EQ( dofNumbers[i] - dofNumbers[i - 1], stride );
  }
}

/**
 * @brief Checks that dof numbers are globally consistent across ranks.
 * @param [in] dofNumbers list of dof numbers, must be sorted
 */
void checkGlobalOrdering( arrayView1d< globalIndex const > const & dofNumbers )
{
  array1d<globalIndex> localDofBounds( 2 );
  localDofBounds[0] = dofNumbers.empty() ? -1 : dofNumbers.front();
  localDofBounds[1] = dofNumbers.empty() ? -1 : dofNumbers.back();

  array1d<globalIndex> globalDofBounds;
  MpiWrapper::allGather( localDofBounds.toViewConst(), globalDofBounds );

  if( MpiWrapper::Comm_rank() == 0 )
  {
    localIndex const newSize =
      std::distance( globalDofBounds.begin(), std::remove( globalDofBounds.begin(), globalDofBounds.end(), -1 ));
    globalDofBounds.resize( newSize );

    // Check that lower ranks have strictly lower dof numbers
    EXPECT_TRUE( std::is_sorted( globalDofBounds.begin(), globalDofBounds.end() ) );
    checkUniqueness( globalDofBounds );
  }
}

/**
 * @brief Test fixture for all non-typed (LAI independent) DofManager tests.
 */
class DofManagerIndicesTest : public DofManagerTestBase
{
protected:

  struct FieldDesc
  {
    string name;
    DofManager::Location location;
    localIndex components;
    std::vector<string> regions = {};
  };

  void testIndices( std::vector<FieldDesc> fields );
};

/**
 * @brief General test function for single-field dof indices.
 * @param location location of the field
 * @param fieldName name of the field
 * @param numComp number of components
 * @param numDofIndicesExpected expected global number of dof indices
 * @param regions list of support regions (empty = whole domain)
 */
void DofManagerIndicesTest::testIndices( std::vector<FieldDesc> fields )
{
  for( FieldDesc & f : fields )
  {
    dofManager.addField( f.name, f.location, f.components, getRegions( mesh, f.regions ) );
  }
  dofManager.reorderByRank();

  array1d<globalIndex> allDofNumbers;
  localIndex lastNumComp = -1;

  for( FieldDesc & f : fields )
  {
    array1d<globalIndex> dofNumbers;
    string const key = dofManager.getKey( f.name );
    string_array const regions = getRegions( mesh, f.regions );
    switch( f.location )
    {
      case DofManager::Location::Elem:
        checkLocalDofNumbers< DofManager::Location::Elem >( mesh, key, regions, dofNumbers );
        break;
      case DofManager::Location::Face:
        checkLocalDofNumbers< DofManager::Location::Face >( mesh, key, regions, dofNumbers );
        break;
      case DofManager::Location::Node:
        checkLocalDofNumbers< DofManager::Location::Node >( mesh, key, regions, dofNumbers );
        break;
      default:
        GEOSX_ERROR( "Unsupported" );
    }
    std::sort( dofNumbers.begin(), dofNumbers.end() );

    if( !dofNumbers.empty() )
    {
      checkUniqueness( dofNumbers );
      checkStride( dofNumbers, f.components );
      if( lastNumComp >= 0 )
      {
        EXPECT_EQ( dofNumbers.front(), allDofNumbers.back() + lastNumComp );
      }
      allDofNumbers.insert( allDofNumbers.size(), dofNumbers.data(), dofNumbers.size() );
    }
    lastNumComp = f.components;
  }
  checkGlobalOrdering( allDofNumbers );
}

/**
 * @brief Check dof indexing for a node-based field with 3 dof/node.
 */
TEST_F( DofManagerIndicesTest, Node_Full )
{
  testIndices( { { "displacement", DofManager::Location::Node, 3 } } );
}

/**
 * @brief Check dof indexing for a node-based field with 3 dof/node with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Node_Partial )
{
  testIndices( { { "displacement", DofManager::Location::Node, 3, { "region1", "region3", "region4" } } } );
}

/**
 * @brief Check dof indexing for a cell-based field with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Elem_Full )
{
  testIndices( { { "pressure", DofManager::Location::Elem, 2 } } );
}

/**
 * @brief Check dof indexing for a cell-based field with 2 dof/cell with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Elem_Partial )
{
  testIndices( { { "pressure", DofManager::Location::Elem, 2, { "region1", "region3", "region4" } } } );
}

/**
 * @brief Check dof indexing for a typical FVM system with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Face_Full )
{
  testIndices( { { "flux", DofManager::Location::Face, 2 } } );
}

/**
 * @brief Check dof indexing for a typical FVM system with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Face_Partial )
{
  testIndices( { { "flux", DofManager::Location::Face, 2, { "region1", "region3", "region4" } } } );
}

/**
 * @brief Check dof indexing for a typical FEM/FVM system.
 */
TEST_F( DofManagerIndicesTest, Node_Elem_Full )
{
  testIndices( { { "displacement", DofManager::Location::Node, 3, {} },
                 { "pressure",     DofManager::Location::Elem, 2, {} } });
}

/**
 * @brief Check dof indexing for a typical FEM/FVM system with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Node_Elem_Partial )
{
  testIndices( { { "displacement", DofManager::Location::Node, 3, { "region1", "region3", "region4" } },
                 { "pressure",     DofManager::Location::Elem, 2, { "region1", "region2", "region4" } } });
}

/**
 * @brief Check dof indexing for a typical Hybrid FVM system.
 */
TEST_F( DofManagerIndicesTest, Face_Elem_Full )
{
  testIndices( { { "flux",     DofManager::Location::Face, 2, {} },
                 { "pressure", DofManager::Location::Elem, 2, {} } });
}

/**
 * @brief Check dof indexing for a typical Hybrid FVM system with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Face_Elem_Partial )
{
  testIndices( { { "flux",     DofManager::Location::Face, 2, { "region1", "region3", "region4" } },
                 { "pressure", DofManager::Location::Elem, 2, { "region1", "region3", "region4" } } });
}

/**
 * @brief Test fixture for all typed (LAI dependent) DofManager tests.
 * @tparam LAI linear algebra interface type
 */
template<typename LAI>
class DofManagerSparsityTest : public DofManagerTestBase
{
protected:

  using PatternFunc = void (*)( MeshLevel const * const mesh,
                                string const & dofIndexKey,
                                string_array const & regions,
                                localIndex const numComp,
                                typename LAI::ParallelMatrix & sparsity );

  using CoupledPatternFunc = void (*)( MeshLevel const * const mesh,
                                       string const & dofIndexKey1,
                                       string const & dofIndexKey2,
                                       string_array const & regions,
                                       localIndex const numComp1,
                                       localIndex const numComp2,
                                       typename LAI::ParallelMatrix & sparsity );

  struct FieldDesc
  {
    string name;
    DofManager::Location location;
    DofManager::Connectivity connectivity;
    localIndex components;
    PatternFunc makePattern;
    std::vector<string> regions = {};
  };

  struct CouplingDesc
  {
    DofManager::Connectivity connectivity;
    CoupledPatternFunc makeCouplingPattern;
    bool symmetric = true;
    std::vector<string> regions = {};
  };

  void testPattern( std::vector<FieldDesc> fields,
                    std::map< std::pair<string, string>, CouplingDesc > couplings = {} );
};

TYPED_TEST_CASE_P( DofManagerSparsityTest );

template<typename LAI>
void DofManagerSparsityTest<LAI>::testPattern( std::vector<FieldDesc> fields,
                                               std::map< std::pair<string, string>, CouplingDesc > couplings )
{
  using Matrix = typename LAI::ParallelMatrix;

  for( FieldDesc const & f : fields )
  {
    dofManager.addField( f.name, f.location, f.components, getRegions( mesh, f.regions ) );
    dofManager.addCoupling( f.name, f.name, f.connectivity );
  }
  for( auto const & entry : couplings )
  {
    std::pair<string, string> const & fieldNames = entry.first;
    CouplingDesc const & c = entry.second;
    dofManager.addCoupling( fieldNames.first, fieldNames.second, c.connectivity, getRegions( mesh, c.regions), c.symmetric );
  }
  dofManager.reorderByRank();

  // Create a sparsity pattern via regular face loop
  localIndex numLocalDof = 0;
  localIndex numCompTotal = 0;
  for( FieldDesc const & f : fields )
  {
    string_array const regions = getRegions( mesh, f.regions );
    localIndex numLocalObj = 0;
    switch( f.location )
    {
      case DofManager::Location::Elem:
        numLocalObj = countLocalObjects< DofManager::Location::Elem >( mesh, regions );
        break;
      case DofManager::Location::Face:
        numLocalObj = countLocalObjects< DofManager::Location::Face >( mesh, regions );
        break;
      case DofManager::Location::Node:
        numLocalObj = countLocalObjects< DofManager::Location::Node >( mesh, regions );
        break;
      default:
        GEOSX_ERROR( "Unsupported" );
    }
    numLocalDof += numLocalObj * f.components;
    numCompTotal += f.components;
  }

  Matrix pattern, patternExpected;

  // Create a sparsity pattern via DofManager
  pattern.createWithLocalSize( dofManager.numLocalDofs(), dofManager.numLocalDofs(), 27 * numCompTotal, MPI_COMM_GEOSX );
  this->dofManager.setSparsityPattern( pattern );

  patternExpected.createWithLocalSize( numLocalDof, numLocalDof, 27 * numCompTotal, MPI_COMM_GEOSX );

  for( FieldDesc const & f : fields )
  {
    f.makePattern( mesh,
                   dofManager.getKey( f.name ),
                   getRegions( mesh, f.regions ),
                   f.components,
                   patternExpected );
  }
  for( auto const & entry : couplings )
  {
    std::pair<string, string> const & names = entry.first;
    CouplingDesc const & c = entry.second;

    FieldDesc const & f1 = *std::find_if( fields.begin(), fields.end(),
                                          [&]( auto const & f ){ return f.name == names.first; } );
    FieldDesc const & f2 = *std::find_if( fields.begin(), fields.end(),
                                          [&]( auto const & f ){ return f.name == names.second; } );

    c.makeCouplingPattern( mesh,
                           dofManager.getKey( f1.name ),
                           dofManager.getKey( f2.name ),
                           getRegions( mesh, c.regions ),
                           f1.components,
                           f2.components,
                           patternExpected );
  }
  patternExpected.close();

  // Compare the sparsity patterns
  pattern.set( 1.0 );
  patternExpected.set( 1.0 );
  compareMatrices( pattern, patternExpected, 0.0, 0.0 );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, TPFA_Full )
{
  TestFixture::testPattern( { { "pressure",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::Face,
                                2, makeSparsityTPFA<TypeParam> } } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, TPFA_Partial )
{
  TestFixture::testPattern( { { "pressure",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::Face,
                                2, makeSparsityTPFA<TypeParam>,
                                { "region1", "region3", "region4" } } } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_Full )
{
  TestFixture::testPattern( { { "displacement",
                                DofManager::Location::Node,
                                DofManager::Connectivity::Elem,
                                3, makeSparsityFEM<TypeParam> } } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_Partial )
{
  TestFixture::testPattern( { { "displacement",
                                DofManager::Location::Node,
                                DofManager::Connectivity::Elem,
                                3, makeSparsityFEM<TypeParam>,
                                { "region1", "region3", "region4" } } } );
}

/**
 * @brief Compare mass matrix sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, Mass_Full )
{
  TestFixture::testPattern( { { "mass",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::None,
                                2, makeSparsityMass<TypeParam> } } );
}

/**
 * @brief Compare mass matrix sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, Mass_Partial )
{
  TestFixture::testPattern( { { "mass",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::None,
                                2, makeSparsityMass<TypeParam>,
                                { "region1", "region3", "region4" } } } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, Flux_Full )
{
  TestFixture::testPattern( { { "flux",
                                DofManager::Location::Face,
                                DofManager::Connectivity::Elem,
                                2, makeSparsityFlux<TypeParam> } } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, Flux_Partial )
{
  TestFixture::testPattern( { { "flux",
                                DofManager::Location::Face,
                                DofManager::Connectivity::Elem,
                                2, makeSparsityFlux<TypeParam>,
                                { "region1", "region3", "region4" } } } );
}

/**
 * @brief Compare a mixed FEM/TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_TPFA_Full )
{
  TestFixture::testPattern( { { "displacement",
                                DofManager::Location::Node,
                                DofManager::Connectivity::Elem,
                                3, makeSparsityFEM<TypeParam> },
                              { "pressure",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::Face,
                                2, makeSparsityTPFA<TypeParam> }
                            },
                            { { { "displacement", "pressure" },
                                { DofManager::Connectivity::Elem,
                                  makeSparsityFEM_FVM<TypeParam>,
                                  true } }
                            } );

}

/**
 * @brief Compare a mixed FEM/TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_TPFA_Partial )
{
  TestFixture::testPattern( { { "displacement",
                                DofManager::Location::Node,
                                DofManager::Connectivity::Elem,
                                3, makeSparsityFEM<TypeParam>,
                                { "region1", "region3", "region4" } },
                              { "pressure",
                                DofManager::Location::Elem,
                                DofManager::Connectivity::Face,
                                2, makeSparsityTPFA<TypeParam>,
                                { "region1", "region2", "region4" }}
                            },
                            { { { "displacement", "pressure" },
                                { DofManager::Connectivity::Elem,
                                  makeSparsityFEM_FVM<TypeParam>,
                                  true,
                                  { "region4" } } }
                            } );

}

REGISTER_TYPED_TEST_CASE_P( DofManagerSparsityTest,
                            TPFA_Full,
                            TPFA_Partial,
                            FEM_Full,
                            FEM_Partial,
                            Mass_Full,
                            Mass_Partial,
                            Flux_Full,
                            Flux_Partial,
                            FEM_TPFA_Full,
                            FEM_TPFA_Partial );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_CASE_P( Trilinos, DofManagerSparsityTest, TrilinosInterface );
#endif

#ifdef GEOSX_USE_PETSC
// Does not work. Remove this comment when fixed.
//INSTANTIATE_TYPED_TEST_CASE_P( Petsc, DofManagerSparsityTest, PetscInterface );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_CASE_P( Hypre, DofManagerSparsityTest, HypreInterface );
#endif

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
