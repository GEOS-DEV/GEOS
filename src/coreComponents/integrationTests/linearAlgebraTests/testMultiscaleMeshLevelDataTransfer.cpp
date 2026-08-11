/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2026 Lawrence Livermore National Security LLC
 * Copyright (c) 2026 GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMultiscaleMeshLevelDataTransfer.cpp
 * @brief Tests transfer of the standard array types through multiscale mesh levels.
 */

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/multiscale/mesh/MeshData.hpp"
#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/DomainPartition.hpp"
#include "integrationTests/linearAlgebraTests/testDofManagerUtils.hpp"

#include <array>
#include <memory>

namespace geos
{
namespace testing
{

namespace
{

char const * const xmlInput =
  R"xml(
  <Problem>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{C3D8}"
                    xCoords="{0, 1}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{2}"
                    ny="{2}"
                    nz="{2}"
                    cellBlockNames="{block}"/>
    </Mesh>
    <ElementRegions>
      <CellElementRegion name="region" cellBlocks="{block}" materialList="{}"/>
    </ElementRegions>
  </Problem>
  )xml";

/**
 * Keep this list synchronized with the layouts in types::StandardArrays. Expanding
 * the operations directly in the test body avoids adding a templated test helper
 * instantiation for every supported array type to coverage's function denominator.
 */
#define GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, VALUE_TYPE ) \
  OP( VALUE_TYPE, 1, RAJA::PERM_I )                         \
  OP( VALUE_TYPE, 2, RAJA::PERM_IJ )                        \
  OP( VALUE_TYPE, 2, RAJA::PERM_JI )                        \
  OP( VALUE_TYPE, 3, RAJA::PERM_IJK )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_IKJ )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_JIK )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_JKI )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_KIJ )                       \
  OP( VALUE_TYPE, 3, RAJA::PERM_KJI )                       \
  OP( VALUE_TYPE, 4, RAJA::PERM_IJKL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IJLK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IKJL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_IKLJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_ILJK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_ILKJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JIKL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JILK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JKIL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JKLI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JLIK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_JLKI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KIJL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KILJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KJIL )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KJLI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KLIJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_KLJI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LIJK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LIKJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LJIK )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LJKI )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LKIJ )                      \
  OP( VALUE_TYPE, 4, RAJA::PERM_LKJI )

#define GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, VALUE_TYPE ) \
  OP( VALUE_TYPE, 1, RAJA::PERM_I )                             \
  OP( VALUE_TYPE, 2, RAJA::PERM_IJ )                            \
  OP( VALUE_TYPE, 2, RAJA::PERM_JI )

#define GEOS_FOR_EACH_STANDARD_ARRAY( OP )               \
  GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, real32 )        \
  GEOS_FOR_EACH_REAL_STANDARD_ARRAY( OP, real64 )        \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, integer )   \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, localIndex ) \
  GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY( OP, globalIndex )

class MultiscaleMeshLevelDataTransferTest : public ::testing::Test
{
protected:

  MultiscaleMeshLevelDataTransferTest()
    : m_state( std::make_unique< CommandLineOptions >() ),
    m_domain( m_state.getProblemManager().getDomainPartition() ),
    m_fineMesh( "fine" ),
    m_coarseMesh( "coarse" )
  {
    setupProblemFromXML( &m_state.getProblemManager(), xmlInput );

    DofManager supportManager( "multiscaleTransferSupport" );
    supportManager.setDomain( m_domain );
    stdVector< DofManager::FieldSupport > const support =
    { { "mesh", "Level0", { "region" } } };
    supportManager.addField( "field",
                             FieldLocation::Elem,
                             1,
                             support );
    m_fineMesh.buildFineMesh( m_domain, supportManager.support( "field" ) );

    LinearSolverParameters::Multiscale::Coarsening coarsening;
    coarsening.partitionType = LinearSolverParameters::Multiscale::Coarsening::PartitionType::cartesian;
    coarsening.ratio.resize( 3 );
    coarsening.ratio[0] = 2;
    coarsening.ratio[1] = 2;
    coarsening.ratio[2] = 2;
    m_coarseMesh.buildCoarseMesh( m_fineMesh,
                                  coarsening,
                                  { "xneg", "xpos", "yneg", "ypos", "zneg", "zpos" } );
  }

  GeosxState m_state;
  DomainPartition & m_domain;
  multiscale::MeshLevel m_fineMesh;
  multiscale::MeshLevel m_coarseMesh;
};

TEST_F( MultiscaleMeshLevelDataTransferTest, standardArraysReachOriginalMesh )
{
  stdVector< string > fieldNames;
  localIndex fieldIndex = 0;

  #define GEOS_REGISTER_TRANSFER_ARRAY( VALUE_TYPE, NDIM, PERMUTATION )                                  \
    do                                                                                                  \
    {                                                                                                   \
      using ArrayType = Array< VALUE_TYPE, NDIM, PERMUTATION >;                                         \
      string const fieldName = "standardArray_" + std::to_string( fieldIndex );                         \
      fieldNames.emplace_back( fieldName );                                                             \
                                                                                                      \
      ArrayType & fineCellField =                                                                        \
        m_fineMesh.cellManager().registerWrapper< ArrayType >( fieldName ).reference();                  \
      std::array< localIndex, NDIM > fineCellDims;                                                       \
      fineCellDims.fill( 2 );                                                                            \
      fineCellDims[0] = m_fineMesh.cellManager().size();                                                \
      fineCellField.resize( NDIM, fineCellDims.data() );                                                \
      for( localIndex i = 0; i < fineCellField.size(); ++i )                                            \
      {                                                                                                 \
        fineCellField.data()[i] = static_cast< VALUE_TYPE >( 256 + 512 * fieldIndex + i );              \
      }                                                                                                 \
                                                                                                      \
      ArrayType & fineNodeField =                                                                        \
        m_fineMesh.nodeManager().registerWrapper< ArrayType >( fieldName ).reference();                  \
      std::array< localIndex, NDIM > fineNodeDims;                                                       \
      fineNodeDims.fill( 2 );                                                                            \
      fineNodeDims[0] = m_fineMesh.nodeManager().size();                                                \
      fineNodeField.resize( NDIM, fineNodeDims.data() );                                                \
      for( localIndex i = 0; i < fineNodeField.size(); ++i )                                            \
      {                                                                                                 \
        fineNodeField.data()[i] = static_cast< VALUE_TYPE >( 32768 + 512 * fieldIndex + i );            \
      }                                                                                                 \
                                                                                                      \
      ArrayType & cellField =                                                                            \
        m_coarseMesh.cellManager().registerWrapper< ArrayType >( fieldName ).reference();                \
      std::array< localIndex, NDIM > cellDims;                                                           \
      cellDims.fill( 2 );                                                                                \
      cellDims[0] = m_coarseMesh.cellManager().size();                                                   \
      cellField.resize( NDIM, cellDims.data() );                                                         \
      for( localIndex i = 0; i < cellField.size(); ++i )                                                 \
      {                                                                                                 \
        cellField.data()[i] = static_cast< VALUE_TYPE >( 1024 + 512 * fieldIndex + i );                  \
      }                                                                                                 \
                                                                                                      \
      ArrayType & nodeField =                                                                            \
        m_coarseMesh.nodeManager().registerWrapper< ArrayType >( fieldName ).reference();                \
      std::array< localIndex, NDIM > nodeDims;                                                           \
      nodeDims.fill( 2 );                                                                                \
      nodeDims[0] = m_coarseMesh.nodeManager().size();                                                   \
      nodeField.resize( NDIM, nodeDims.data() );                                                         \
      for( localIndex i = 0; i < nodeField.size(); ++i )                                                 \
      {                                                                                                 \
        nodeField.data()[i] = static_cast< VALUE_TYPE >( 65536 + 512 * fieldIndex + i );                 \
      }                                                                                                 \
      ++fieldIndex;                                                                                     \
    } while( false );

  GEOS_FOR_EACH_STANDARD_ARRAY( GEOS_REGISTER_TRANSFER_ARRAY )
#undef GEOS_REGISTER_TRANSFER_ARRAY

  ASSERT_EQ( fieldIndex, 75 );
  m_fineMesh.writeCellData( fieldNames );
  m_fineMesh.writeNodeData( fieldNames );
  m_coarseMesh.writeCellData( fieldNames );
  m_coarseMesh.writeNodeData( fieldNames );

  geos::MeshLevel & sourceMesh = m_domain.getMeshBody( "mesh" ).getMeshLevel( "Level0" );
  ElementSubRegionBase & sourceSubRegion = sourceMesh.getElemManager().getRegion( 0 ).getSubRegion( 0 );

  arrayView1d< localIndex const > const coarseCellIndex =
    m_fineMesh.cellManager().getField< fields::multiscale::CoarseCellLocalIndex >();
  arrayView1d< localIndex const > const origRegion =
    m_fineMesh.cellManager().getField< fields::multiscale::OrigElementRegion >();
  arrayView1d< localIndex const > const origSubRegion =
    m_fineMesh.cellManager().getField< fields::multiscale::OrigElementSubRegion >();
  arrayView1d< localIndex const > const origCellIndex =
    m_fineMesh.cellManager().getField< fields::multiscale::OrigElementIndex >();
  for( localIndex i = 0; i < m_fineMesh.cellManager().size(); ++i )
  {
    ASSERT_EQ( origRegion[i], 0 );
    ASSERT_EQ( origSubRegion[i], 0 );
  }

  arrayView1d< localIndex const > const fineNodeIndex =
    m_coarseMesh.nodeManager().getField< fields::multiscale::FineNodeLocalIndex >();
  arrayView1d< localIndex const > const origNodeIndex =
    m_fineMesh.nodeManager().getField< fields::multiscale::OrigNodeIndex >();

  fieldIndex = 0;
  #define GEOS_VERIFY_TRANSFERRED_ARRAY( VALUE_TYPE, NDIM, PERMUTATION )                                 \
    do                                                                                                  \
    {                                                                                                   \
      using ArrayType = Array< VALUE_TYPE, NDIM, PERMUTATION >;                                         \
      string const fieldName = "standardArray_" + std::to_string( fieldIndex );                         \
      string const fineOutputName = m_fineMesh.name() + '_' + fieldName;                                \
      string const coarseOutputName = m_coarseMesh.name() + '_' + fieldName;                            \
                                                                                                      \
      ArrayType const & fineCellField =                                                                 \
        m_fineMesh.cellManager().getReference< ArrayType >( fieldName );                                \
      ArrayType const & fineOutputCellField =                                                           \
        sourceSubRegion.getReference< ArrayType >( fineOutputName );                                    \
      for( integer dim = 1; dim < NDIM; ++dim )                                                         \
      {                                                                                                 \
        EXPECT_EQ( fineOutputCellField.size( dim ), fineCellField.size( dim ) );                         \
      }                                                                                                 \
      for( localIndex i = 0; i < m_fineMesh.cellManager().size(); ++i )                                 \
      {                                                                                                 \
        LvArray::forValuesInSliceWithIndices( fineCellField[i],                                         \
                                              [&]( auto const & expectedValue, auto const ... indices )  \
{                                                                                               \
  EXPECT_EQ( fineOutputCellField( origCellIndex[i], indices ... ), expectedValue );              \
} );                                                                                            \
      }                                                                                                 \
                                                                                                      \
      ArrayType const & coarseCellField =                                                               \
        m_coarseMesh.cellManager().getReference< ArrayType >( fieldName );                              \
      ArrayType const & coarseOutputCellField =                                                         \
        sourceSubRegion.getReference< ArrayType >( coarseOutputName );                                  \
      for( integer dim = 1; dim < NDIM; ++dim )                                                         \
      {                                                                                                 \
        EXPECT_EQ( coarseOutputCellField.size( dim ), coarseCellField.size( dim ) );                     \
      }                                                                                                 \
      for( localIndex i = 0; i < m_fineMesh.cellManager().size(); ++i )                                 \
      {                                                                                                 \
        LvArray::forValuesInSliceWithIndices( coarseCellField[coarseCellIndex[i]],                       \
                                              [&]( auto const & expectedValue, auto const ... indices )  \
{                                                                                               \
  EXPECT_EQ( coarseOutputCellField( origCellIndex[i], indices ... ), expectedValue );            \
} );                                                                                            \
      }                                                                                                 \
                                                                                                      \
      ArrayType const & fineNodeField =                                                                 \
        m_fineMesh.nodeManager().getReference< ArrayType >( fieldName );                                \
      ArrayType const & fineOutputNodeField =                                                           \
        sourceMesh.getNodeManager().getReference< ArrayType >( fineOutputName );                        \
      for( integer dim = 1; dim < NDIM; ++dim )                                                         \
      {                                                                                                 \
        EXPECT_EQ( fineOutputNodeField.size( dim ), fineNodeField.size( dim ) );                         \
      }                                                                                                 \
      for( localIndex i = 0; i < m_fineMesh.nodeManager().size(); ++i )                                 \
      {                                                                                                 \
        LvArray::forValuesInSliceWithIndices( fineNodeField[i],                                         \
                                              [&]( auto const & expectedValue, auto const ... indices )  \
{                                                                                               \
  EXPECT_EQ( fineOutputNodeField( origNodeIndex[i], indices ... ), expectedValue );              \
} );                                                                                            \
      }                                                                                                 \
                                                                                                      \
      ArrayType const & coarseNodeField =                                                               \
        m_coarseMesh.nodeManager().getReference< ArrayType >( fieldName );                              \
      ArrayType const & coarseOutputNodeField =                                                         \
        sourceMesh.getNodeManager().getReference< ArrayType >( coarseOutputName );                      \
      for( integer dim = 1; dim < NDIM; ++dim )                                                         \
      {                                                                                                 \
        EXPECT_EQ( coarseOutputNodeField.size( dim ), coarseNodeField.size( dim ) );                     \
      }                                                                                                 \
      for( localIndex i = 0; i < m_coarseMesh.nodeManager().size(); ++i )                               \
      {                                                                                                 \
        LvArray::forValuesInSliceWithIndices( coarseNodeField[i],                                       \
                                              [&]( auto const & expectedValue, auto const ... indices )  \
{                                                                                               \
  EXPECT_EQ( coarseOutputNodeField( origNodeIndex[fineNodeIndex[i]], indices ... ),              \
             expectedValue );                                                                   \
} );                                                                                            \
      }                                                                                                 \
      ++fieldIndex;                                                                                     \
    } while( false );

  GEOS_FOR_EACH_STANDARD_ARRAY( GEOS_VERIFY_TRANSFERRED_ARRAY )
#undef GEOS_VERIFY_TRANSFERRED_ARRAY

  EXPECT_EQ( fieldIndex, 75 );
}

#undef GEOS_FOR_EACH_STANDARD_ARRAY
#undef GEOS_FOR_EACH_INTEGRAL_STANDARD_ARRAY
#undef GEOS_FOR_EACH_REAL_STANDARD_ARRAY

} // namespace
} // namespace testing
} // namespace geos

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
