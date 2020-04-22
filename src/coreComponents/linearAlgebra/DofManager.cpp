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
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationOps.hpp"
#include "mesh/MeshLevel.hpp"

#include "DofManagerHelpers.hpp"

#include <numeric>

namespace geosx
{

using namespace dataRepository;

DofManager::DofManager( string name )
  : m_name( std::move( name ) ),
  m_domain( nullptr ),
  m_mesh( nullptr ),
  m_reordered( false )
{
  initializeDataStructure();
}

void DofManager::initializeDataStructure()
{
  // we pre-allocate an oversized array to store connectivity type
  // instead of resizing it dynamically as fields are added.
  m_coupling.resize( MAX_FIELDS, MAX_FIELDS );
}

void DofManager::clear()
{
  // deallocate index arrays from the mesh
  for( FieldDescription const & field : m_fields )
  {
    removeIndexArray( field );
  }

  // delete internal data
  m_fields.clear();
  m_coupling.clear();

  initializeDataStructure();

  m_reordered = false;
}

void DofManager::setMesh( DomainPartition * const domain,
                          localIndex const meshLevelIndex,
                          localIndex const meshBodyIndex )
{
  // TODO: this should be m_domain != domain
  if( m_domain != nullptr )
  {
    // Domain is changed! Delete old data structure and create new
    clear();
  }
  m_domain = domain;
  m_mesh = m_domain->getMeshBody( meshBodyIndex )->getMeshLevel( meshLevelIndex );
}

localIndex DofManager::getFieldIndex( string const & name ) const
{
  GEOSX_ASSERT_MSG( fieldExists( name ), "DofManager: field does not exist: " << name );
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  return std::distance( m_fields.begin(), it );
}

bool DofManager::fieldExists( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  return it != m_fields.end();
}

string DofManager::getKey( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].key;
}

globalIndex DofManager::numGlobalDofs( string const & fieldName ) const
{
  if( !fieldName.empty() )
  {
    return m_fields[getFieldIndex( fieldName )].numGlobalDof;
  }
  else
  {
    return std::accumulate( m_fields.begin(), m_fields.end(), 0,
                            []( globalIndex const & n, FieldDescription const & f ) { return n + f.numGlobalDof; } );
  }
}

localIndex DofManager::numLocalDofs( string const & fieldName ) const
{
  if( !fieldName.empty() )
  {
    return m_fields[getFieldIndex( fieldName )].numLocalDof;
  }
  else
  {
    return std::accumulate( m_fields.begin(), m_fields.end(), 0,
                            []( globalIndex const & n, FieldDescription const & f ) { return n + f.numLocalDof; } );
  }
}

localIndex DofManager::rankOffset( string const & fieldName ) const
{
  if( !fieldName.empty() )
  {
    return m_fields[getFieldIndex( fieldName )].rankOffset;
  }
  else
  {
    return std::accumulate( m_fields.begin(), m_fields.end(), 0,
                            []( globalIndex const & n, FieldDescription const & f ) { return n + f.rankOffset; } );
  }
}

localIndex DofManager::numComponents( string const & fieldName ) const
{
  if( !fieldName.empty() )
  {
    return m_fields[getFieldIndex( fieldName )].numComponents;
  }
  else
  {
    return std::accumulate( m_fields.begin(), m_fields.end(), 0,
                            []( globalIndex const & n, FieldDescription const & f ) { return n + f.numComponents; } );
  }
}

void DofManager::createIndexArray( FieldDescription & field )
{
  bool const success =
    LocationSwitch( field.location, [&]( auto const loc )
  {
    Location constexpr LOC = decltype(loc)::value;
    using helper = IndexArrayHelper< globalIndex, LOC >;

    // 0. register index arrays
    helper::template create<>( m_mesh, field.key, field.docstring, field.regions );
    typename helper::Accessor & indexArray = helper::get( m_mesh, field.key );

    // step 1. loop over all active regions, determine number of local mesh objects
    localIndex numLocalNodes = 0;
    forMeshLocation< LOC, false >( m_mesh, field.regions, [&]( auto const locIdx )
    {
      helper::reference( indexArray, locIdx ) = field.numComponents * numLocalNodes++;
    } );
    field.numLocalDof = field.numComponents * numLocalNodes;

    // step 2. gather row counts across ranks
    field.rankOffset = MpiWrapper::PrefixSum< globalIndex >( field.numLocalDof );

    field.numGlobalDof = field.rankOffset + field.numLocalDof;
    MpiWrapper::Broadcast( field.numGlobalDof, MpiWrapper::Comm_size() - 1 );

    // step 3. adjust local dof offsets to reflect processor offset
    forMeshLocation< LOC, false >( m_mesh, field.regions, [&]( auto const locIdx )
    {
      helper::reference( indexArray, locIdx ) += field.rankOffset;
    } );

    // step 4. synchronize across ranks
    std::map< string, string_array > fieldNames;
    fieldNames[ MeshHelper< LOC >::syncObjName ].push_back( field.key );

    CommunicationTools::
      SynchronizeFields( fieldNames, m_mesh,
                         m_domain->getNeighbors() );
  } );
  GEOSX_ERROR_IF( !success, "Invalid location type: " << static_cast< int >(field.location) );
}

void DofManager::removeIndexArray( FieldDescription const & field )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    Location constexpr LOC = decltype(loc)::value;
    IndexArrayHelper< globalIndex, LOC >::template remove<>( m_mesh, field.key, field.regions );
  } );
}

// Just an interface to allow only three parameters
void DofManager::addField( string const & fieldName,
                           Location const location )
{
  addField( fieldName, location, 1, array1d< string >() );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addField( string const & fieldName,
                           Location const location,
                           localIndex const components )
{
  addField( fieldName, location, components, array1d< string >() );
}

// Just another interface to allow four parameters (no components)
void DofManager::addField( string const & fieldName,
                           Location const location,
                           arrayView1d< string const > const & regions )
{
  addField( fieldName, location, 1, regions );
}

// The real function, allowing the creation of self-connected blocks
void DofManager::addField( string const & fieldName,
                           Location const location,
                           localIndex const components,
                           arrayView1d< string const > const & regions )
{
  GEOSX_ERROR_IF( m_reordered, "Cannot add fields after reorderByRank() has been called." );
  GEOSX_ERROR_IF( fieldExists( fieldName ), "Requested field name '" << fieldName << "' already exists." );
  GEOSX_ERROR_IF( m_fields.size() >= MAX_FIELDS, "Limit on DofManager's MAX_NUM_FIELDS exceeded." );

  localIndex fieldIndex = m_fields.size();
  m_fields.resize( fieldIndex + 1 );

  string suffix;
  for( string const & regionName : regions )
  {
    suffix.append( "_" + regionName );
  }

  FieldDescription & field = m_fields.back();

  field.name = fieldName;
  field.location = location;
  field.regions.resize( regions.size() );
  std::copy( regions.begin(), regions.end(), field.regions.begin() );
  field.numComponents = components;
  field.key = m_name + '_' + fieldName + "_dofIndex" + suffix;
  field.docstring = fieldName + " DoF indices";

  if( components > 1 )
  {
    field.docstring += " (with " + std::to_string( components ) + "-component blocks)";
  }

  ElementRegionManager * const elemManager = m_mesh->getElemManager();
  if( field.regions.empty() )
  {
    elemManager->forElementRegions( [&]( ElementRegionBase const & region )
    {
      field.regions.push_back( region.getName() );
    } );
  }
  else
  {
    for( string const & regionName : field.regions )
    {
      GEOSX_ERROR_IF( elemManager->GetRegion( regionName ) == nullptr, "Element region not found: " << regionName );
    }
  }

  // sort and remove duplicates
  std::sort( field.regions.begin(), field.regions.end() );
  auto const end_it = std::unique( field.regions.begin(), field.regions.end() );
  field.regions.resize( std::distance( field.regions.begin(), end_it ) );

  m_coupling[fieldIndex][fieldIndex].regions = field.regions;

  // based on location, allocate an index array for this field
  createIndexArray( field );

  // determine field's global offset
  if( fieldIndex > 0 )
  {
    FieldDescription & prev = m_fields[fieldIndex - 1];
    field.blockOffset = prev.blockOffset + prev.numGlobalDof;
  }
  else
  {
    field.blockOffset = 0;
  }
  field.globalOffset = -1; // to be calculated in reorderByRank()
}

namespace
{

/**
 * @brief Implements a mesh loop that populates a local connector-to-location pattern.
 * @tparam LOC type of mesh locations
 * @tparam CONN type of mesh connector
 * @tparam SUBREGIONTYPES types of subregions to loop over
 *
 * The algorithm used requires that a connector-to-location map is present and populated in
 * the mesh data structure. We loop over connectors first, then visit adjacent locations from
 * each connector and thus populate sparsity pattern row-by-row. The numbering of connectors
 * is implicit, i.e. implied by the looping order, thus we don't have to compute/store it.
 *
 * The sparsity pattern thus constructed should only be contracted with another sparsity pattern that
 * was built with the same type of connectors and list of region names and types, otherwise connector
 * numbering will be incompatible between the two.
 */
template< DofManager::Location LOC, DofManager::Location CONN, typename ... SUBREGIONTYPES >
struct ConnLocPatternBuilder
{
  static void build( MeshLevel * const mesh,
                     DofManager::FieldDescription const & field,
                     array1d< string > const & regions,
                     localIndex const rowOffset,
                     LvArray::SparsityPattern< globalIndex > & connLocPattern )
  {
    using ArrayHelper = IndexArrayHelper< globalIndex const, LOC >;
    typename ArrayHelper::Accessor dofIndexArray = ArrayHelper::get( mesh, field.key );
    auto connIdxPrev = MeshHelper< CONN >::invalid_local_index;

    localIndex connectorCount = -1;

    forMeshLocation< CONN, LOC, true, SUBREGIONTYPES... >( mesh, regions,
                                                           [&]( auto const & connIdx,
                                                                auto const & locIdx,
                                                                localIndex const GEOSX_UNUSED_PARAM( k ) )
    {
      if( connIdx != connIdxPrev )
      {
        ++connectorCount;
        connIdxPrev = connIdx;
      }

      globalIndex const dofOffset = ArrayHelper::value( dofIndexArray, locIdx );
      if( dofOffset >= 0 )
      {
        for( localIndex c = 0; c < field.numComponents; ++c )
        {
          connLocPattern.insertNonZero( rowOffset + connectorCount, dofOffset + c );
        }
      }
    } );
  }
};

/**
 * @brief A specialization of ConnLocPatternBuilder for elements connected through edges.
 * @tparam SUBREGIONTYPES types of subregions to loop over
 *
 * This is required because edge-to-element map does not exist in the mesh. Therefore, we have
 * to invert the loop order and go through elements first. This requires us to create a unique
 * connector numbering for edges in the regions of interest first - which comes with a memory
 * and runtime overhead, so we only use this method when absolutely have to, i.e. in this case.
 */
template< typename ... SUBREGIONTYPES >
struct ConnLocPatternBuilder< DofManager::Location::Elem, DofManager::Location::Edge, SUBREGIONTYPES... >
{
  static void build( MeshLevel * const mesh,
                     DofManager::FieldDescription const & field,
                     array1d< string > const & regions,
                     localIndex const rowOffset,
                     LvArray::SparsityPattern< globalIndex > & connLocPattern )
  {
    DofManager::Location constexpr ELEM  = DofManager::Location::Elem;
    DofManager::Location constexpr EDGE = DofManager::Location::Edge;

    EdgeManager const * const edgeManager = mesh->getEdgeManager();
    ElementRegionManager const * const elemManager = mesh->getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofIndex =
      elemManager->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( field.key );

    array1d< localIndex > edgeConnectorIndex( edgeManager->size() );
    edgeConnectorIndex = -1;

    localIndex edgeCount = 0;
    forMeshLocation< EDGE, true, SUBREGIONTYPES... >( mesh, regions, [&] ( localIndex const edgeIdx )
    {
      edgeConnectorIndex[edgeIdx] = edgeCount++;
    } );

    forMeshLocation< ELEM, EDGE, true, SUBREGIONTYPES... >( mesh, regions,
                                                            [&]( auto const & elemIdx,
                                                                 auto const & edgeIdx,
                                                                 localIndex const GEOSX_UNUSED_PARAM( k ) )
    {
      globalIndex const dofOffset = dofIndex[elemIdx[0]][elemIdx[1]][elemIdx[2]];
      if( dofOffset >= 0 )
      {
        for( localIndex c = 0; c < field.numComponents; ++c )
        {
          connLocPattern.insertNonZero( rowOffset + edgeConnectorIndex[edgeIdx], dofOffset + c );
        }
      }
    } );
  }
};

template< DofManager::Location LOC, DofManager::Location CONN >
void makeConnLocPattern( MeshLevel * const mesh,
                         DofManager::FieldDescription const & field,
                         array1d< string > const & regions,
                         LvArray::SparsityPattern< globalIndex > & connLocPattern )
{
  using Loc = DofManager::Location;

  // 1. Count number of connectors in order to size the sparsity pattern
  localIndex const numConnectors = countMeshObjects< CONN, true >( mesh, regions );
  localIndex numConnectorsTotal = numConnectors;

  // TPFA+fracture special case
  if( LOC == Loc::Elem && CONN == Loc::Face )
  {
    numConnectorsTotal += countMeshObjects< Loc::Edge, true, FaceElementSubRegion >( mesh, regions );
  }

  // 2. Resize the pattern appropriately
  localIndex const numEntriesPerRow = MeshIncidence< CONN, LOC >::max * field.numComponents;
  connLocPattern.resize( numConnectorsTotal,
                         field.numGlobalDof,
                         numEntriesPerRow < field.numGlobalDof ? numEntriesPerRow : field.numGlobalDof );

  // 3. Populate the local CL pattern
  ConnLocPatternBuilder< LOC, CONN >::build( mesh, field, regions, 0, connLocPattern );

  // TPFA+fracture special case
  if( LOC == Loc::Elem && CONN == Loc::Face )
  {
    ConnLocPatternBuilder< Loc::Elem, Loc::Edge, FaceElementSubRegion >::build( mesh, field, regions, numConnectors, connLocPattern );
  }
}

} // namespace

template< typename MATRIX >
void DofManager::setSparsityPatternFromStencil( MATRIX & pattern,
                                                localIndex const fieldIndex ) const
{
  FieldDescription const & field = m_fields[fieldIndex];
  CouplingDescription const & coupling = m_coupling[fieldIndex][fieldIndex];
  localIndex const NC = field.numComponents;

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumber =
    m_mesh->getElemManager()->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( field.key );

  array1d< globalIndex > rowIndices( NC );
  array1d< globalIndex > colIndices( NC );
  array2d< real64 > values( NC, NC );
  values = 1.0;

  // 1. Insert diagonal blocks, in case there are elements not included in stencil
  // (e.g. a single fracture element not connected to any other)
  forMeshLocation< Location::Elem, false >( m_mesh, field.regions,
                                            [&]( auto const & elemIdx )
  {
    for( localIndex c = 0; c < NC; ++c )
    {
      rowIndices[c] = dofNumber[elemIdx[0]][elemIdx[1]][elemIdx[2]] + c;
    }
    pattern.insert( rowIndices, rowIndices, values );
  } );

  // 2. Assemble diagonal and off-diagonal blocks for elements in stencil
  MATRIX * const pattern_ptr = &pattern;
  coupling.stencils->forAllStencils( [&]( auto const & stencil )
  {
    using StenciType = typename std::decay< decltype( stencil ) >::type;
    constexpr localIndex maxNumFluxElems = StenciType::NUM_POINT_IN_FLUX;
    constexpr localIndex maxStencilSize = StenciType::MAX_STENCIL_SIZE;

    typename StenciType::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    typename StenciType::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    typename StenciType::IndexContainerViewConstType const & sei = stencil.getElementIndices();

    rowIndices.reserve( maxNumFluxElems * NC );
    colIndices.reserve( maxStencilSize * NC );
    values.reserve( maxNumFluxElems * NC * maxStencilSize * NC );

    forAll< serialPolicy >( stencil.size(), [&]( localIndex iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      localIndex const stencilSize  = numFluxElems;

      rowIndices.resize( numFluxElems * NC );
      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        for( localIndex c = 0; c < NC; ++c )
        {
          rowIndices[i * NC + c] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + c;
        }
      }

      colIndices.resize( stencilSize * NC );
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        for( localIndex c = 0; c < NC; ++c )
        {
          colIndices[i * NC + c] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + c;
        }
      }

      values.resize( numFluxElems * NC, stencilSize * NC );
      values = 1.0;

      pattern_ptr->insert( rowIndices, colIndices, values );
    } );
  } );
}

template< typename MATRIX >
void DofManager::setSparsityPatternOneBlock( MATRIX & pattern,
                                             localIndex const rowFieldIndex,
                                             localIndex const colFieldIndex ) const
{
  GEOSX_ASSERT( rowFieldIndex >= 0 );
  GEOSX_ASSERT( colFieldIndex >= 0 );

  FieldDescription const & rowField = m_fields[rowFieldIndex];
  FieldDescription const & colField = m_fields[colFieldIndex];

  Connector conn = m_coupling[rowFieldIndex][colFieldIndex].connector;
  if( rowFieldIndex == colFieldIndex && conn == Connector::Stencil )
  {
    setSparsityPatternFromStencil( pattern, rowFieldIndex );
    return;
  }

  LvArray::SparsityPattern< globalIndex > connLocRow( 0, 0, 0 ), connLocCol( 0, 0, 0 );

  localIndex maxDofRow = 0;
  LocationSwitch( rowField.location, static_cast< Location >( conn ),
                  [&]( auto const locType, auto const connType )
  {
    Location constexpr LOC  = decltype(locType)::value;
    Location constexpr CONN = decltype(connType)::value;

    makeConnLocPattern< LOC, CONN >( m_mesh,
                                     rowField,
                                     m_coupling[rowFieldIndex][colFieldIndex].regions,
                                     connLocRow );

    maxDofRow = MeshIncidence< CONN, LOC >::max * rowField.numComponents;
  } );

  localIndex maxDofCol = 0;
  if( colFieldIndex == rowFieldIndex )
  {
    connLocCol = connLocRow; // TODO avoid copying
    maxDofCol = maxDofRow;
  }
  else
  {
    LocationSwitch( colField.location, static_cast< Location >( conn ),
                    [&]( auto const locType, auto const connType )
    {
      Location constexpr LOC = decltype(locType)::value;
      Location constexpr CONN = decltype(connType)::value;

      makeConnLocPattern< LOC, CONN >( m_mesh,
                                       colField,
                                       m_coupling[rowFieldIndex][colFieldIndex].regions,
                                       connLocCol );

      maxDofCol = MeshIncidence< CONN, LOC >::max * colField.numComponents;
    } );
  }
  GEOSX_ASSERT_EQ( connLocRow.numRows(), connLocCol.numRows() );

  array1d< globalIndex > dofIndicesRow( maxDofRow );
  array1d< globalIndex > dofIndicesCol( maxDofCol );
  array2d< real64 > values( maxDofRow, maxDofCol );

  // Perform assembly/multiply patterns
  for( localIndex irow = 0; irow < connLocRow.numRows(); ++irow )
  {
    localIndex const numDofRow = connLocRow.numNonZeros( irow );
    dofIndicesRow.resize( numDofRow );
    for( localIndex j = 0; j < numDofRow; ++j )
    {
      dofIndicesRow[j] = connLocRow.getColumns( irow )[j];
    }

    localIndex const numDofCol = connLocCol.numNonZeros( irow );
    dofIndicesCol.resize( numDofCol );
    for( localIndex j = 0; j < numDofCol; ++j )
    {
      dofIndicesCol[j] = connLocCol.getColumns( irow )[j];
    }

    values.resize( numDofRow, numDofCol );
    values = 1.0;
    pattern.insert( dofIndicesRow, dofIndicesCol, values );
  }
}

// Create the sparsity pattern (location-location). Low level interface
template< typename MATRIX >
void DofManager::setSparsityPattern( MATRIX & matrix,
                                     bool const closePattern ) const
{
  GEOSX_ERROR_IF( !m_reordered, "Cannot set monolithic sparsity pattern before reorderByRank() has been called." );

  matrix.open();
  for( localIndex blockRow = 0; blockRow < m_fields.size(); ++blockRow )
  {
    for( localIndex blockCol = 0; blockCol < m_fields.size(); ++blockCol )
    {
      setSparsityPatternOneBlock( matrix, blockRow, blockCol );
    }
  }
  if( closePattern )
  {
    matrix.close();
  }
}

// Create the sparsity pattern (location-location). High level interface
template< typename MATRIX >
void DofManager::setSparsityPattern( MATRIX & matrix,
                                     string const & rowFieldName,
                                     string const & colFieldName,
                                     bool const closePattern ) const
{
  GEOSX_ERROR_IF( m_reordered, "Cannot set single block sparsity pattern after reorderByRank() has been called." );
  setSparsityPatternOneBlock( matrix, getFieldIndex( rowFieldName ), getFieldIndex( colFieldName ) );
  if( closePattern )
  {
    matrix.close();
  }
}

namespace
{

template< typename FIELD_OP, typename POLICY, typename VECTOR >
void vectorToFieldImpl( VECTOR const & vector,
                        ObjectManagerBase & manager,
                        string const & dofKey,
                        string const & fieldName,
                        real64 const scalingFactor,
                        localIndex const dofOffset,
                        localIndex const loComp,
                        localIndex const hiComp )
{
  arrayView1d< globalIndex const > const & indexArray = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = manager.ghostRank();

  real64 const * localVector = vector.extractLocalVector();

  WrapperBase * const wrapper = manager.getWrapperBase( fieldName );
  GEOSX_ASSERT( wrapper != nullptr );

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( std::type_index( wrapper->get_typeid() ) ),
                                  false,
                                  [&]( auto arrayInstance,
                                       auto GEOSX_UNUSED_PARAM( dataTypeInstance ) )
  {
    using ArrayType = decltype( arrayInstance );
    Wrapper< ArrayType > & view = Wrapper< ArrayType >::cast( *wrapper );
    traits::ViewType< ArrayType > field = view.reference();

    forAll< POLICY >( indexArray.size(), [=] ( localIndex const i )
    {
      if( ghostRank[i] < 0 && indexArray[i] >= 0 )
      {
        localIndex const lid = indexArray[i] - dofOffset;
        GEOSX_ASSERT( lid >= 0 );
        for( localIndex c = loComp; c < hiComp; ++c )
        {
          FIELD_OP::template SpecifyFieldValue( field,
                                                i,
                                                integer_conversion< integer >( c - loComp ),
                                                scalingFactor * localVector[lid + c] );
        }
      }
    } );
  } );
}

template< typename FIELD_OP, typename POLICY, typename VECTOR >
void fieldToVectorImpl( VECTOR & vector,
                        ObjectManagerBase const & manager,
                        string const & dofKey,
                        string const & fieldName,
                        real64 const GEOSX_UNUSED_PARAM( scalingFactor ),
                        localIndex const dofOffset,
                        localIndex const loComp,
                        localIndex const hiComp )
{
  arrayView1d< globalIndex const > const & indexArray = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = manager.ghostRank();

  real64 * localVector = vector.extractLocalVector();

  WrapperBase const * const wrapper = manager.getWrapperBase( fieldName );
  GEOSX_ASSERT( wrapper != nullptr );

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( std::type_index( wrapper->get_typeid() ) ),
                                  false,
                                  [&]( auto arrayInstance,
                                       auto GEOSX_UNUSED_PARAM( dataTypeInstance ) )
  {
    using ArrayType = decltype( arrayInstance );
    Wrapper< ArrayType > const & view = Wrapper< ArrayType >::cast( *wrapper );
    traits::ViewTypeConst< ArrayType > field = view.reference();

    forAll< POLICY >( indexArray.size(), [=] ( localIndex const i )
    {
      if( ghostRank[i] < 0 && indexArray[i] >= 0 )
      {
        localIndex const lid = indexArray[i] - dofOffset;
        GEOSX_ASSERT( lid >= 0 );
        for( localIndex c = loComp; c < hiComp; ++c )
        {
          FIELD_OP::template ReadFieldValue( field,
                                             i,
                                             integer_conversion< integer >( c - loComp ),
                                             localVector[lid + c] );
        }
      }
    } );
  } );
}

}

template< typename FIELD_OP, typename POLICY, typename VECTOR >
void DofManager::vectorToField( VECTOR const & vector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                real64 const scalingFactor,
                                localIndex const loCompIndex,
                                localIndex const hiCompIndex ) const
{
  FieldDescription const & fieldDesc = m_fields[getFieldIndex( srcFieldName )];

  localIndex const loComp = loCompIndex;
  localIndex const hiComp = (hiCompIndex >= 0) ? hiCompIndex : fieldDesc.numComponents;
  GEOSX_ASSERT( loComp >= 0 && hiComp <= fieldDesc.numComponents && loComp < hiComp );

  if( fieldDesc.location == Location::Elem )
  {
    m_mesh->getElemManager()->forElementSubRegions< ElementSubRegionBase >( fieldDesc.regions,
                                                                            [&]( localIndex const,
                                                                                 ElementSubRegionBase & subRegion )
    {
      vectorToFieldImpl< FIELD_OP, POLICY >( vector,
                                             subRegion,
                                             fieldDesc.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             loComp,
                                             hiComp );
    } );
  }
  else
  {
    vectorToFieldImpl< FIELD_OP, POLICY >( vector,
                                           getObjectManager( fieldDesc.location, m_mesh ),
                                           fieldDesc.key,
                                           dstFieldName,
                                           scalingFactor,
                                           rankOffset(),
                                           loComp,
                                           hiComp );
  }
}

// Copy values from DOFs to nodes
template< typename VECTOR >
void DofManager::copyVectorToField( VECTOR const & vector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    localIndex const loCompIndex,
                                    localIndex const hiCompIndex ) const
{
  vectorToField< FieldSpecificationEqual, parallelHostPolicy >( vector,
                                                                srcFieldName,
                                                                dstFieldName,
                                                                scalingFactor,
                                                                loCompIndex,
                                                                hiCompIndex );
}

// Copy values from DOFs to nodes
template< typename VECTOR >
void DofManager::addVectorToField( VECTOR const & vector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   real64 const scalingFactor,
                                   localIndex const loCompIndex,
                                   localIndex const hiCompIndex ) const
{
  vectorToField< FieldSpecificationAdd, parallelHostPolicy >( vector,
                                                              srcFieldName,
                                                              dstFieldName,
                                                              scalingFactor,
                                                              loCompIndex,
                                                              hiCompIndex );
}

template< typename FIELD_OP, typename POLICY, typename VECTOR >
void DofManager::fieldToVector( VECTOR & vector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                real64 const scalingFactor,
                                localIndex const loCompIndex,
                                localIndex const hiCompIndex ) const
{
  FieldDescription const & fieldDesc = m_fields[getFieldIndex( srcFieldName )];

  localIndex const loComp = loCompIndex;
  localIndex const hiComp = (hiCompIndex >= 0) ? hiCompIndex : fieldDesc.numComponents;
  GEOSX_ASSERT( loComp >= 0 && hiComp <= fieldDesc.numComponents && loComp < hiComp );

  if( fieldDesc.location == Location::Elem )
  {
    m_mesh->getElemManager()->forElementSubRegions< ElementSubRegionBase >( fieldDesc.regions,
                                                                            [&]( localIndex const,
                                                                                 ElementSubRegionBase const & subRegion )
    {
      fieldToVectorImpl< FIELD_OP, POLICY >( vector,
                                             subRegion,
                                             fieldDesc.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             loComp,
                                             hiComp );
    } );
  }
  else
  {
    fieldToVectorImpl< FIELD_OP, POLICY >( vector,
                                           getObjectManager( fieldDesc.location, m_mesh ),
                                           fieldDesc.key,
                                           dstFieldName,
                                           scalingFactor,
                                           rankOffset(),
                                           loComp,
                                           hiComp );
  }
}

// Copy values from nodes to DOFs
template< typename VECTOR >
void DofManager::copyFieldToVector( VECTOR & vector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    localIndex const loCompIndex,
                                    localIndex const hiCompIndex ) const
{
  fieldToVector< FieldSpecificationEqual, parallelHostPolicy >( vector,
                                                                srcFieldName,
                                                                dstFieldName,
                                                                scalingFactor,
                                                                loCompIndex,
                                                                hiCompIndex );
}

// Copy values from nodes to DOFs
template< typename VECTOR >
void DofManager::addFieldToVector( VECTOR & vector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   real64 const scalingFactor,
                                   localIndex const loCompIndex,
                                   localIndex const hiCompIndex ) const
{
  fieldToVector< FieldSpecificationAdd, parallelHostPolicy >( vector,
                                                              srcFieldName,
                                                              dstFieldName,
                                                              scalingFactor,
                                                              loCompIndex,
                                                              hiCompIndex );
}

// Just an interface to allow only three parameters
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity )
{
  string_array dummy_regions;
  addCoupling( rowFieldName, colFieldName, connectivity, dummy_regions, true );
}

// Just another interface to allow four parameters (no symmetry)
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity,
                              arrayView1d< string const > const & regions )
{
  addCoupling( rowFieldName, colFieldName, connectivity, regions, true );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity,
                              bool const symmetric )
{
  string_array dummy_regions;
  addCoupling( rowFieldName, colFieldName, connectivity, dummy_regions, symmetric );
}

// The real function, allowing the creation of coupling blocks
void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity,
                              arrayView1d< string const > const & regions,
                              bool const symmetric )
{
  localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
  localIndex const colFieldIndex = getFieldIndex( colFieldName );

  // Check if already defined
  if( m_coupling[rowFieldIndex][colFieldIndex].connector != Connector::None )
  {
    GEOSX_ERROR( "addCoupling: coupling already defined with another connectivity" );
    return;
  }

  // get row/col field regions
  string_array const & rowRegions = m_fields[rowFieldIndex].regions;
  string_array const & colRegions = m_fields[colFieldIndex].regions;
  string_array & regionList = m_coupling[rowFieldIndex][colFieldIndex].regions;

  if( regions.empty() )
  {
    // Populate with regions common between two fields
    regionList.resize( std::min( rowRegions.size(), colRegions.size() ) );
    string_array::iterator it = std::set_intersection( rowRegions.begin(), rowRegions.end(),
                                                       colRegions.begin(), colRegions.end(),
                                                       regionList.begin() );
    regionList.resize( std::distance( regionList.begin(), it ) );
  }
  else
  {
    // Sort alphabetically and remove possible duplicates in user input
    regionList.resize( regions.size() );
    std::copy( regions.begin(), regions.end(), regionList.begin() );
    std::sort( regionList.begin(), regionList.end() );
    regionList.resize( std::distance( regionList.begin(), std::unique( regionList.begin(), regionList.end() ) ) );

    // Check that both fields exist on all regions in the list
    for( string const & regionName : regionList )
    {
      GEOSX_ERROR_IF( std::find( rowRegions.begin(), rowRegions.end(), regionName ) == rowRegions.end(),
                      "Region '" << regionName << "' does not belong to the domain of field '" << rowFieldName << "'" );
      GEOSX_ERROR_IF( std::find( colRegions.begin(), colRegions.end(), regionName ) == colRegions.end(),
                      "Region '" << regionName << "' does not belong to the domain of field '" << colFieldName << "'" );
    }
  }

  // save field's connectivity type (rowField to colField)
  m_coupling[rowFieldIndex][colFieldIndex].connector = connectivity;

  if( connectivity == Connector::None && rowFieldIndex == colFieldIndex )
  {
    m_coupling[rowFieldIndex][colFieldIndex].connector = static_cast< Connector >( m_fields[rowFieldIndex].location );
  }

  // Set connectivity with active symmetry flag
  if( symmetric && colFieldIndex != rowFieldIndex )
  {
    m_coupling[colFieldIndex][rowFieldIndex].connector = connectivity;
    m_coupling[colFieldIndex][rowFieldIndex].regions = regionList;
  }
}

void DofManager::addCoupling( string const & fieldName,
                              FluxApproximationBase const * stencils )
{
  localIndex const fieldIndex = getFieldIndex( fieldName );
  FieldDescription const & field = m_fields[fieldIndex];

  GEOSX_ERROR_IF( field.location != Location::Elem, "Field must be supported on elements in order to use stencil sparsity" );
  GEOSX_ERROR_IF( stencils == nullptr, "Must pass a valid pointer to a flux approximation object" );

  CouplingDescription & coupling = m_coupling[fieldIndex][fieldIndex];
  coupling.connector = Connector::Stencil;
  coupling.regions = field.regions;
  coupling.stencils = stencils;
}

void DofManager::reorderByRank()
{
  GEOSX_LAI_ASSERT( !m_reordered );

  // update field offsets to account for renumbering
  globalIndex dofOffset = rankOffset();
  for( FieldDescription & field : m_fields )
  {
    field.globalOffset = dofOffset;
    dofOffset += field.numLocalDof;
  }

  std::map< string, string_array > fieldToSync;

  // adjust index arrays for owned locations
  for( FieldDescription const & field : m_fields )
  {
    globalIndex const adjustment = field.globalOffset - field.rankOffset;

    LocationSwitch( field.location, [&]( auto const loc )
    {
      Location constexpr LOC = decltype(loc)::value;
      using ArrayHelper = IndexArrayHelper< globalIndex, LOC >;

      typename ArrayHelper::Accessor indexArray = ArrayHelper::get( m_mesh, field.key );

      forMeshLocation< LOC, false >( m_mesh, field.regions, [&]( auto const locIdx )
      {
        ArrayHelper::reference( indexArray, locIdx ) += adjustment;
      } );

      fieldToSync[MeshHelper< LOC >::syncObjName].push_back( field.key );
    } );
  }

  // synchronize index arrays for all fields across ranks
  CommunicationTools::
    SynchronizeFields( fieldToSync, m_mesh,
                       m_domain->getNeighbors() );

  m_reordered = true;
}

template< typename MATRIX >
void DofManager::makeRestrictor( string const & fieldName,
                                 MATRIX & restrictor,
                                 MPI_Comm const comm,
                                 bool const transpose,
                                 localIndex const loCompIndex,
                                 localIndex const hiCompIndex ) const
{
  GEOSX_ERROR_IF( !m_reordered, "Cannot make restrictors before reorderByRank() has been called." );

  FieldDescription const & field = m_fields[getFieldIndex( fieldName )];

  localIndex const loComp = loCompIndex;
  localIndex const hiComp = (hiCompIndex >= 0) ? hiCompIndex : field.numComponents;
  GEOSX_ASSERT( loComp >= 0 && hiComp <= field.numComponents && loComp < hiComp );

  localIndex const numComp = hiComp - loComp;
  localIndex const numLoc = field.numLocalDof / field.numComponents;

  globalIndex rowOffset = field.rankOffset / field.numComponents * numComp - loComp;
  globalIndex colOffset = field.globalOffset;
  localIndex rowStride = numComp;
  localIndex colStride = field.numComponents;
  localIndex rowSize = field.numLocalDof / field.numComponents * numComp;
  localIndex colSize = numLocalDofs();

  if( transpose )
  {
    std::swap( rowOffset, colOffset );
    std::swap( rowStride, colStride );
    std::swap( rowSize, colSize );
  }

  restrictor.createWithLocalSize( rowSize, colSize, 1, comm );
  restrictor.open();
  for( localIndex i = 0; i < numLoc; ++i )
  {
    for( localIndex c = loComp; c < hiComp; ++c )
    {
      restrictor.insert( rowOffset + i * rowStride + c, colOffset + i * colStride + c, 1.0 );
    }
  }
  restrictor.close();
}

// Print the coupling table on screen
void DofManager::printFieldInfo( std::ostream & os ) const
{
  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) == 0 )
  {
    localIndex numFields = m_fields.size();

    os << "Fields:" << std::endl;
    os << " # | " << std::setw( 20 ) << "name" << " | " << "comp" << " | " << "N global DOF" << std::endl;
    os << "---+----------------------+------+-------------" << std::endl;
    for( localIndex i = 0; i < numFields; ++i )
    {
      FieldDescription const & f = m_fields[i];
      os << ' ' << i << " | "
         << std::setw( 20 ) << f.name << " | "
         << std::setw( 4 ) << f.numComponents << " | "
         << std::setw( 12 ) << f.numGlobalDof << std::endl;
    }
    os << "---+----------------------+------+-------------" << std::endl;

    os << std::endl << "Connectivity:" << std::endl;
    for( localIndex i = 0; i < numFields; ++i )
    {
      for( localIndex j = 0; j < numFields; ++j )
      {
        switch( m_coupling[i][j].connector )
        {
          case Connector::Elem:
          {
            os << " E ";
          }
          break;
          case Connector::Face:
          {
            os << " F ";
          }
          break;
          case Connector::Edge:
          {
            os << " G ";
          }
          break;
          case Connector::Node:
          {
            os << " N ";
          }
          break;
          case Connector::None:
          {
            os << "   ";
          }
          break;
          case Connector::Stencil:
          {
            os << " S ";
          }
          break;
          default:
            GEOSX_ERROR( "Invalid connector type: " << static_cast< int >( m_coupling[i][j].connector ) );
        }
        if( j < numFields - 1 )
        {
          os << "|";
        }
      }
      os << std::endl;
      if( i < numFields - 1 )
      {
        for( localIndex j = 0; j < numFields - 1; ++j )
        {
          os << "---|";
        }
        os << "---";
      }
      os << std::endl;
    }
  }
}

#define MAKE_DOFMANAGER_METHOD_INST( LAI ) \
  template void DofManager::setSparsityPattern( LAI::ParallelMatrix &, \
                                                bool const ) const; \
  template void DofManager::setSparsityPattern( LAI::ParallelMatrix &, \
                                                string const &, \
                                                string const &, \
                                                bool const ) const; \
  template void DofManager::copyVectorToField( LAI::ParallelVector const &, \
                                               string const &, \
                                               string const &, \
                                               real64 const, \
                                               localIndex const, \
                                               localIndex const ) const; \
  template void DofManager::addVectorToField( LAI::ParallelVector const &, \
                                              string const &, \
                                              string const &, \
                                              real64 const, \
                                              localIndex const, \
                                              localIndex const ) const; \
  template void DofManager::copyFieldToVector( LAI::ParallelVector &, \
                                               string const &, \
                                               string const &, \
                                               real64 const, \
                                               localIndex const, \
                                               localIndex const ) const; \
  template void DofManager::addFieldToVector( LAI::ParallelVector &, \
                                              string const &, \
                                              string const &, \
                                              real64 const, \
                                              localIndex const, \
                                              localIndex const ) const; \
  template void DofManager::makeRestrictor( string const &, \
                                            LAI::ParallelMatrix &, \
                                            MPI_Comm const, \
                                            bool const, \
                                            localIndex const, \
                                            localIndex const ) const;

#ifdef GEOSX_USE_TRILINOS
MAKE_DOFMANAGER_METHOD_INST( TrilinosInterface )
#endif

#ifdef GEOSX_USE_HYPRE
MAKE_DOFMANAGER_METHOD_INST( HypreInterface )
#endif

#ifdef GEOSX_USE_PETSC
MAKE_DOFMANAGER_METHOD_INST( PetscInterface )
#endif

} // namespace geosx
