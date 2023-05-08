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
 * @file Mesh.cpp
 */

#include "MeshLevel.hpp"

#include "common/TypeDispatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/multiscale/mesh/coarsening/Coarsening.hpp"
#include "MeshData.hpp"
#include "MeshUtils.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

using namespace dataRepository;

namespace multiscale
{

MeshLevel::MeshLevel( string const & name )
  : m_name( name ),
  m_root( m_name, m_rootNode ),
  m_cellManager( m_name + "_cellManager", &m_root ),
  m_nodeManager( m_name + "_nodeManager", &m_root )
{}

namespace
{

localIndex assignNodeLocalIndices( NodeManager & nodeManager,
                                   string const & localIndexKey,
                                   string const & globalIndexKey,
                                   globalIndex const rankOffset,
                                   localIndex const numLocalNodes )
{
  arrayView1d< globalIndex const > const nodeGlobalIndex = nodeManager.getReference< array1d< globalIndex > >( globalIndexKey ).toViewConst();
  arrayView1d< localIndex > const nodeLocalIndex = nodeManager.getReference< array1d< localIndex > >( localIndexKey ).toView();
  localIndex numPresentNodes = numLocalNodes;

  // Assign new local indices s.t. locally owned objects come first, followed by ghosted ones.
  // We do this by looping and filtering over all nodes of the original mesh, as this is the primary use case.
  // It might pessimize the case of a small flow region in a larger mechanics mesh, but it's less critical.
  // Note: there doesn't seem to be a good way to do this in parallel ... or does it?
  forAll< serialPolicy >( nodeManager.size(), [&]( localIndex const k )
  {
    globalIndex const newGlobalIndex = nodeGlobalIndex[k];
    GEOS_ERROR_IF_LT( newGlobalIndex, 0 );
    if( newGlobalIndex >= 0 )
    {
      localIndex const localIndexTmp = static_cast< localIndex >( newGlobalIndex - rankOffset );
      localIndex const newLocalIndex = ( 0 <= localIndexTmp && localIndexTmp < numLocalNodes ) ? localIndexTmp : numPresentNodes++;
      nodeLocalIndex[k] = newLocalIndex;
    }
  } );

  return numPresentNodes;
}

void copyNodeData( NodeManager const & nodeManager,
                   string const & localIndexKey,
                   string const & globalIndexKey,
                   MeshObjectManager & msNodeManager )
{
  arrayView1d< localIndex const > const nodeLocalIndex = nodeManager.getReference< array1d< localIndex > >( localIndexKey );

  meshUtils::fillArrayByDstIndex< parallelHostPolicy >( nodeLocalIndex,
                                                        msNodeManager.getField< fields::multiscale::OrigNodeIndex >().toView(),
                                                        []( auto _ ){ return _; } );
  meshUtils::fillArrayByDstIndex< parallelHostPolicy >( nodeLocalIndex,
                                                        msNodeManager.localToGlobalMap().toView(),
                                                        nodeManager.getReference< array1d< globalIndex > >( globalIndexKey ) );
  meshUtils::fillArrayByDstIndex< parallelHostPolicy >( nodeLocalIndex,
                                                        msNodeManager.ghostRank().toView(),
                                                        nodeManager.ghostRank().toViewConst() );
  meshUtils::fillArrayByDstIndex< parallelHostPolicy >( nodeLocalIndex,
                                                        msNodeManager.isExternal().toView(),
                                                        nodeManager.isExternal().toViewConst() );
  meshUtils::fillArrayByDstIndex< parallelHostPolicy >( nodeLocalIndex,
                                                        msNodeManager.getDomainBoundaryIndicator().toView(),
                                                        nodeManager.getDomainBoundaryIndicator().toViewConst() );
  msNodeManager.constructGlobalToLocalMap();
  msNodeManager.setMaxGlobalIndex();
}

void populateNodeManager( string const & localIndexKey,
                          DomainPartition & domain,
                          DofManager::FieldSupport const & support,
                          MeshObjectManager & msNodeManager )
{
  // Use a temporary DofManager to generate new contiguous global index numbering
  string const dofName = "globalIndex";
  DofManager dofManager( "nodeIndexMgr" );
  dofManager.setDomain( domain );
  dofManager.addField( dofName, FieldLocation::Node, 1, { support } );
  dofManager.reorderByRank();

  string const globalIndexKey = dofManager.getKey( dofName );
  geos::MeshLevel & mesh = domain.getMeshBody( support.meshBodyName ).getMeshLevel( support.meshLevelName );
  NodeManager & nodeManager = mesh.getNodeManager();

  localIndex const numLocalNodes = dofManager.numLocalDofs( dofName );
  localIndex const numPresentNodes = assignNodeLocalIndices( nodeManager,
                                                             localIndexKey,
                                                             globalIndexKey,
                                                             dofManager.rankOffset(),
                                                             numLocalNodes );

  // Setup neighbor data
  std::vector< int > const neighborRanks = domain.getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    msNodeManager.addNeighbor( rank );
  }

  // Populate the new node manager, now that we know its total size
  msNodeManager.resize( numPresentNodes );
  msNodeManager.setNumOwnedObjects( numLocalNodes );
  copyNodeData( nodeManager, localIndexKey, globalIndexKey, msNodeManager );
  meshUtils::copySets( nodeManager, localIndexKey, msNodeManager );
  meshUtils::copyNeighborData( nodeManager, localIndexKey, neighborRanks, msNodeManager, meshUtils::mapIndexArray< localIndex > );

  // Remove index fields on the GEOSX mesh
  dofManager.clear();
}

localIndex assignCellLocalIndices( ElementRegionManager & elemManager,
                                   std::set< string > const & regions,
                                   string const & localIndexKey,
                                   string const & globalIndexKey,
                                   globalIndex const rankOffset,
                                   localIndex const numLocalCells )
{
  localIndex numPresentCells = numLocalCells;

  // Assign new local indices s.t. locally owned objects come first, followed by ghosted ones.
  elemManager.forElementSubRegions( regions, [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    // Allocate an element-based field on the GEOSX mesh that will keep a mapping to new cell local indices
    arrayView1d< localIndex > const cellLocalIndex = subRegion.getReference< array1d< localIndex > >( localIndexKey ).toView();
    arrayView1d< globalIndex const > const cellGlobalIndex = subRegion.getReference< array1d< globalIndex > >( globalIndexKey ).toViewConst();

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const ei )
    {
      localIndex const localIndexTmp = static_cast< localIndex >( cellGlobalIndex[ei] - rankOffset );
      cellLocalIndex[ei] = ( 0 <= localIndexTmp && localIndexTmp < numLocalCells ) ? localIndexTmp : numPresentCells++;
    } );
  } );

  return numPresentCells;
}

void copyCellData( ElementRegionManager const & elemManager,
                   std::set< string > const & regions,
                   string const & localIndexKey,
                   string const & globalIndexKey,
                   MeshObjectManager & msCellManager )
{
  elemManager.forElementSubRegionsComplete( regions, [&]( localIndex const,
                                                          localIndex const er,
                                                          localIndex const esr,
                                                          ElementRegionBase const &,
                                                          ElementSubRegionBase const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex = subRegion.getReference< array1d< localIndex > >( localIndexKey ).toViewConst();

    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.getField< fields::multiscale::OrigElementRegion >().toView(),
                                                          [er]( auto ){ return er; } );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.getField< fields::multiscale::OrigElementSubRegion >().toView(),
                                                          [esr]( auto ){ return esr; } );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.getField< fields::multiscale::OrigElementIndex >().toView(),
                                                          []( auto _ ){ return _; } );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.localToGlobalMap().toView(),
                                                          subRegion.getReference< array1d< globalIndex > >( globalIndexKey ) );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.ghostRank().toView(),
                                                          subRegion.ghostRank().toViewConst() );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.isExternal().toView(),
                                                          subRegion.isExternal().toViewConst() );
    meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                          msCellManager.getDomainBoundaryIndicator().toView(),
                                                          subRegion.getDomainBoundaryIndicator().toViewConst() );

    if( subRegion.hasField< fields::StructuredIndex >() )
    {
      arrayView2d< integer const > const origStructIndex = subRegion.getField< fields::StructuredIndex >();
      arrayView2d< integer > const msStructIndex = [&]()
      {
        if( msCellManager.hasWrapper( fields::StructuredIndex::key() ) )
        {
          return msCellManager.getField< fields::StructuredIndex >().toView();
        }
        else
        {
          array2d< integer > & msCartIndexArray =
            msCellManager.registerField< fields::StructuredIndex >( "multiscale::MeshLevel" ).reference();
          msCartIndexArray.resizeDimension< 1 >( origStructIndex.size( 1 ) );
          msCartIndexArray.setValues< parallelHostPolicy >( -1 );
          return msCartIndexArray.toView();
        }
      }();

      meshUtils::fillArrayByDstIndex< parallelHostPolicy >( cellLocalIndex,
                                                            msStructIndex.toView(),
                                                            origStructIndex );
    }
  } );
  msCellManager.constructGlobalToLocalMap();
  msCellManager.setMaxGlobalIndex();
}

void populateCellManager( string const & localIndexKey,
                          DomainPartition & domain,
                          DofManager::FieldSupport const & support,
                          MeshObjectManager & msCellManager )
{
  // Use a temporary DofManager to generate new contiguous global index numbering
  string const dofName = "globalIndex";
  DofManager dofManager( "cellIndexMgr" );
  dofManager.setDomain( domain );
  dofManager.addField( dofName, FieldLocation::Elem, 1, { support } );
  dofManager.reorderByRank();

  string const globalIndexKey = dofManager.getKey( dofName );
  geos::MeshLevel & mesh = domain.getMeshBody( support.meshBodyName ).getMeshLevel( support.meshLevelName );
  ElementRegionManager & elemManager = mesh.getElemManager();

  localIndex const numLocalCells = dofManager.numLocalDofs( dofName );
  localIndex const numPresentCells = assignCellLocalIndices( elemManager,
                                                             support.regionNames,
                                                             localIndexKey,
                                                             globalIndexKey,
                                                             dofManager.rankOffset(),
                                                             numLocalCells );

  // Setup neighbor data
  std::vector< int > const neighborRanks = domain.getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    msCellManager.addNeighbor( rank );
  }

  // Populate the new cell manager, now that we know its total size
  msCellManager.resize( numPresentCells );
  msCellManager.setNumOwnedObjects( numLocalCells );
  copyCellData( elemManager, support.regionNames, localIndexKey, globalIndexKey, msCellManager );
  elemManager.forElementSubRegions( support.regionNames, [&]( localIndex, ElementSubRegionBase const & subRegion )
  {
    meshUtils::copySets( subRegion, localIndexKey, msCellManager );
    meshUtils::copyNeighborData( subRegion, localIndexKey, neighborRanks, msCellManager, meshUtils::mapIndexArray< localIndex > );
  } );

  // Make sure to remove temporary global index field on the GEOSX mesh
  dofManager.clear();
}

void buildNodeToCellMap( string const & cellLocalIndexKey,
                         geos::MeshLevel const & mesh,
                         MeshObjectManager & msNodeManager )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const origNodeIndex = msNodeManager.getField< fields::multiscale::OrigNodeIndex >().toViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex const > > cellLocalIndex =
    elemManager.constructArrayViewAccessor< localIndex, 1 >( cellLocalIndexKey );

  ArrayOfArraysView< localIndex const > const elemRegion = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const elemSubRegion = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const elemIndex = nodeManager.elementList().toViewConst();

  // Count the length of each sub-array
  array1d< localIndex > rowCounts( msNodeManager.size() );
  forAll< parallelHostPolicy >( msNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const k )
  {
    rowCounts[k] = elemIndex.sizeOfArray( origNodeIndex[k] );
  } );

  // Resize
  MeshObjectManager::MapType & nodeToCell = msNodeManager.toDualRelation();
  nodeToCell.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  forAll< parallelHostPolicy >( msNodeManager.size(), [=, nodeToCell = nodeToCell.toView()]( localIndex const k )
  {
    localIndex const nodeIdx = origNodeIndex[k];
    arraySlice1d< localIndex const > const er = elemRegion[nodeIdx];
    arraySlice1d< localIndex const > const es = elemSubRegion[nodeIdx];
    arraySlice1d< localIndex const > const ei = elemIndex[nodeIdx];
    for( localIndex a = 0; a < ei.size(); ++a )
    {
      nodeToCell.insertIntoSet( k, cellLocalIndex[er[a]][es[a]][ei[a]] );
    }
  } );
}

void buildCellToNodeMap( string const & cellLocalIndexKey,
                         string const & nodeLocalIndexKey,
                         geos::MeshLevel const & mesh,
                         DofManager::FieldSupport const & support,
                         MeshObjectManager & msCellManager )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const nodeLocalIndex =
    nodeManager.getReference< array1d< localIndex > >( nodeLocalIndexKey ).toViewConst();

  // Count the length of each row
  array1d< localIndex > rowCounts( msCellManager.size() );
  elemManager.forElementSubRegions( support.regionNames, [&]( localIndex, auto const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex =
      subRegion.template getReference< array1d< localIndex > >( cellLocalIndexKey ).toViewConst();
    auto const elemToNode = subRegion.nodeList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [cellLocalIndex, elemToNode,
                                                     rowCounts = rowCounts.toView()]( localIndex const ei )
    {
      rowCounts[cellLocalIndex[ei]] = elemToNode[ei].size();
    } );
  } );

  // Resize the map
  MeshObjectManager::MapType & cellToNode = msCellManager.toDualRelation();
  cellToNode.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  elemManager.forElementSubRegions( support.regionNames, [&]( localIndex, auto const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex =
      subRegion.template getReference< array1d< localIndex > >( cellLocalIndexKey ).toViewConst();
    auto const elemToNode = subRegion.nodeList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [cellLocalIndex, elemToNode, nodeLocalIndex,
                                                     cellToNode = cellToNode.toView()]( localIndex const ei )
    {
      auto const nodes = elemToNode[ei];
      localIndex const cellIdx = cellLocalIndex[ei];
      for( localIndex a = 0; a < nodes.size(); ++a )
      {
        cellToNode.insertIntoSet( cellIdx, nodeLocalIndex[nodes[a]] );
      }
    } );
  } );
}

} // namespace

void MeshLevel::buildFineMesh( DomainPartition & domain,
                               Span< geos::DofManager::FieldSupport const > const support )
{
  GEOS_MARK_FUNCTION;

  GEOS_ERROR_IF_NE_MSG( support.size(), 1,
                        GEOS_FMT( "[MsRSB] {}: multiple domain bodies/levels not supported", m_name ) );
  m_domain = &domain;
  m_support = support[0];
  geos::MeshLevel & sourceMesh = domain.getMeshBody( m_support.meshBodyName ).getMeshLevel( m_support.meshLevelName );

  string const cellLocalIndexKey = m_name + "_cell_localIndex";
  string const nodeLocalIndexKey = m_name + "_node_localIndex";

  // Allocate fields on multiscale sourceMesh to keep mapping to original elements/nodes
  m_cellManager.registerField< fields::multiscale::OrigElementRegion >( m_name );
  m_cellManager.registerField< fields::multiscale::OrigElementSubRegion >( m_name );
  m_cellManager.registerField< fields::multiscale::OrigElementIndex >( m_name );
  m_nodeManager.registerField< fields::multiscale::OrigNodeIndex >( m_name );

  // Allocate fields on the GEOSX sourceMesh that will keep mappings to new local indices
  sourceMesh.getNodeManager().registerWrapper< array1d< localIndex > >( nodeLocalIndexKey ).setApplyDefaultValue( -1 );
  sourceMesh.getElemManager().forElementSubRegions( m_support.regionNames,
                                                    [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    subRegion.registerWrapper< array1d< localIndex > >( cellLocalIndexKey ).setApplyDefaultValue( -1 );
  } );

  // Generate new contiguous local/global numberings for the target subset of nodes/cells
  populateNodeManager( nodeLocalIndexKey, *m_domain, m_support, m_nodeManager );
  populateCellManager( cellLocalIndexKey, *m_domain, m_support, m_cellManager );

  // Extract appropriately renumbered cell-node maps into new managers
  buildNodeToCellMap( cellLocalIndexKey, sourceMesh, m_nodeManager );
  buildCellToNodeMap( cellLocalIndexKey, nodeLocalIndexKey, sourceMesh, m_support, m_cellManager );

  // Remove local index fields used during construction from GEOSX sourceMesh
  sourceMesh.getNodeManager().deregisterWrapper( nodeLocalIndexKey );
  sourceMesh.getElemManager().forElementSubRegions( m_support.regionNames,
                                                    [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    subRegion.deregisterWrapper( cellLocalIndexKey );
  } );
}

void MeshLevel::buildCoarseMesh( multiscale::MeshLevel & fineMesh,
                                 LinearSolverParameters::Multiscale::Coarsening const & coarse_params,
                                 array1d< string > const & boundaryNodeSets )
{
  GEOS_MARK_FUNCTION;
  m_fineMesh = &fineMesh;
  m_domain = fineMesh.m_domain;
  coarsening::buildCoarseMesh( fineMesh, *this, coarse_params, boundaryNodeSets );
}

void MeshLevel::writeCellData( std::vector< string > const & fieldNames, int depth ) const
{
  if( m_fineMesh )
  {
    writeCellDataCoarse( fieldNames, depth );
  }
  else
  {
    writeCellDataFine( fieldNames, depth );
  }
}

void MeshLevel::writeNodeData( std::vector< string > const & fieldNames, int depth ) const
{
  if( m_fineMesh )
  {
    writeNodeDataCoarse( fieldNames, depth );
  }
  else
  {
    writeNodeDataFine( fieldNames, depth );
  }
}

void MeshLevel::writeCellDataFine( std::vector< string > const & fieldNames, int depth ) const
{
  geos::MeshLevel & sourceMesh = m_domain->getMeshBody( m_support.meshBodyName ).getMeshLevel( m_support.meshLevelName );
  arrayView1d< localIndex const > const origRegion    = m_cellManager.getField< fields::multiscale::OrigElementRegion >();
  arrayView1d< localIndex const > const origSubRegion = m_cellManager.getField< fields::multiscale::OrigElementSubRegion >();
  arrayView1d< localIndex const > const origIndex     = m_cellManager.getField< fields::multiscale::OrigElementIndex >();

  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_cellManager.getWrapperBase( fieldName );
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;
      using ArrayViewType = typename ArrayType::ParentClass;

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const sourceFieldName = ( depth == 0 ? m_name + '_' : "" ) + wrapper.getName();
      sourceMesh.getElemManager().forElementSubRegions( m_support.regionNames,
                                                        [&]( localIndex, ElementSubRegionBase & subRegion )
      {
        subRegion.registerWrapper< ArrayType >( sourceFieldName ).copyWrapperAttributes( wrapper );
        auto & dstField = subRegion.getReference< ArrayType >( sourceFieldName );

        // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
        auto dims = srcField.dimsArray();
        dims[0] = dstField.size( 0 );
        dstField.resize( ArrayType::NDIM, dims.data );
      } );

      auto accessor = sourceMesh.getElemManager().constructViewAccessor< ArrayType, ArrayViewType >( sourceFieldName );
      forAll< parallelHostPolicy >( m_cellManager.size(), [=, dstField = accessor.toNestedView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField[origRegion[ic]][origSubRegion[ic]]( origIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    }, wrapper );
  }
}

void MeshLevel::writeCellDataCoarse( std::vector< string > const & fieldNames, int depth ) const
{
  GEOS_ASSERT( m_fineMesh != nullptr );
  arrayView1d< localIndex const > const coarseCellIndex =
    m_fineMesh->cellManager().getField< fields::multiscale::CoarseCellLocalIndex >();

  std::vector< string > fineFieldNames;
  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_cellManager.getWrapperBase( fieldName );
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const fineFieldName = ( depth == 0 ? m_name + '_' : "" ) + wrapper.getName();
      fineFieldNames.push_back( fineFieldName );
      m_fineMesh->cellManager().registerWrapper< ArrayType >( fineFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = m_fineMesh->cellManager().getReference< ArrayType >( fineFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( dstField.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( dstField[ ic ],
                                              [&]( auto & dstVal, auto const ... indices )
        {
          dstVal = srcField( coarseCellIndex[ic], indices ... );
        } );
      } );
    }, wrapper );
  }

  // Recursively call on finer levels
  m_fineMesh->writeCellData( fineFieldNames, depth + 1 );
}

void MeshLevel::writeNodeDataFine( std::vector< string > const & fieldNames, int depth ) const
{
  geos::MeshLevel & sourceMesh = m_domain->getMeshBody( m_support.meshBodyName ).getMeshLevel( m_support.meshLevelName );
  arrayView1d< localIndex const > const origNodeIndex = m_nodeManager.getField< fields::multiscale::OrigNodeIndex >();

  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_nodeManager.getWrapperBase( fieldName );
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const sourceFieldName = ( depth == 0 ? m_name + '_' : "" ) + wrapper.getName();
      sourceMesh.getNodeManager().registerWrapper< ArrayType >( sourceFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = sourceMesh.getNodeManager().getReference< ArrayType >( sourceFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( m_nodeManager.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField( origNodeIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    }, wrapper );
  }
}

void MeshLevel::writeNodeDataCoarse( std::vector< string > const & fieldNames, int depth ) const
{
  GEOS_ASSERT( m_fineMesh != nullptr );
  arrayView1d< localIndex const > const fineNodeIndex = m_nodeManager.getField< fields::multiscale::FineNodeLocalIndex >();

  std::vector< string > fineFieldNames;
  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_nodeManager.getWrapperBase( fieldName );
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const fineFieldName = ( depth == 0 ? m_name + '_' : "" ) + wrapper.getName();
      fineFieldNames.push_back( fineFieldName );
      m_fineMesh->nodeManager().registerWrapper< ArrayType >( fineFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = m_fineMesh->nodeManager().getReference< ArrayType >( fineFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( m_nodeManager.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField( fineNodeIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    }, wrapper );
  }

  // Recursively call on finer levels
  m_fineMesh->writeNodeData( fineFieldNames, depth + 1 );
}

} // namespace multiscale
} // namespace geos
