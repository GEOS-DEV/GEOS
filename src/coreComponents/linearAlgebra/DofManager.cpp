/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

#include "common/FieldSpecificationOps.hpp"
#include "common/TypeDispatch.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"

#include "DofManagerHelpers.hpp"

#include <numeric>
#include <functional>

namespace geosx
{

using namespace dataRepository;

DofManager::DofManager( string name )
  : m_name( std::move( name ) )
{}

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
  m_reordered = false;
}

void DofManager::setDomain( DomainPartition & domain )
{
  clear();
  m_domain = &domain;
}

localIndex DofManager::getFieldIndex( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  GEOSX_ASSERT_MSG( it != m_fields.end(), "DofManager: field does not exist: " << name );
  return std::distance( m_fields.begin(), it );
}

bool DofManager::fieldExists( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  return it != m_fields.end();
}

string const & DofManager::getKey( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].key;
}

globalIndex DofManager::numGlobalDofs( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numGlobalDof;
}

globalIndex DofManager::numGlobalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numGlobalDof; } );
}

localIndex DofManager::numLocalDofs( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numLocalDof;
}

localIndex DofManager::numLocalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), localIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numLocalDof; } );
}

globalIndex DofManager::rankOffset( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].rankOffset;
}

globalIndex DofManager::rankOffset() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.rankOffset; } );
}

integer DofManager::numComponents( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numComponents;
}

integer DofManager::numComponents() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), integer( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numComponents; } );
}

array1d< integer > DofManager::numComponentsPerField() const
{
  array1d< integer > ret( m_fields.size() );
  std::transform( m_fields.begin(), m_fields.end(), ret.begin(), std::mem_fn( &FieldDescription::numComponents ) );
  return ret;
}

DofManager::Location DofManager::location( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].location;
}

globalIndex DofManager::globalOffset( string const & fieldName ) const
{
  GEOSX_ASSERT_MSG( m_reordered, "Global offset not available until after reorderByRank() has been called." );
  return m_fields[getFieldIndex( fieldName )].globalOffset;
}

namespace
{

template< typename FUNC >
void forMeshSupport( std::vector< DofManager::Regions > const & support,
                     DomainPartition & domain,
                     FUNC && func )
{
  for( DofManager::Regions const & regions : support )
  {
    MeshBody & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );
    func( meshBody, meshLevel, regions.regionNames );
  }
}

template< typename FUNC >
void forMeshSupport( std::vector< DofManager::Regions > const & support,
                     DomainPartition const & domain,
                     FUNC && func )
{
  for( DofManager::Regions const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );
    func( meshBody, meshLevel, regions.regionNames );
  }
}

} // namespace

void DofManager::createIndexArray( FieldDescription const & field )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    localIndex index = 0;

    forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
    {
      Location constexpr LOC = decltype(loc)::value;
      using helper = ArrayHelper< globalIndex, LOC >;

      std::map< string, string_array > fieldNames;
      fieldNames[MeshHelper< LOC >::syncObjName].emplace_back( field.key );

      // register index array
      helper::template create<>( mesh, field.key, field.docstring, regions );
      typename helper::Accessor indexArray = helper::get( mesh, field.key );

      // populate index array using a sequential counter
      forMeshLocation< LOC, false, serialPolicy >( mesh, regions, [&]( auto const locIdx )
      {
        helper::reference( indexArray, locIdx ) = field.rankOffset + field.numComponents * index++;
      } );

      // synchronize across ranks
      CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, m_domain->getNeighbors(), false );
    } );
  } );
}

void DofManager::removeIndexArray( FieldDescription const & field )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
    {
      Location constexpr LOC = decltype(loc)::value;
      ArrayHelper< globalIndex, LOC >::template remove<>( mesh, field.key, regions );
    } );
  } );
}

void DofManager::computeFieldDimensions( localIndex const fieldIndex )
{
  FieldDescription & field = m_fields[fieldIndex];

  // determine number of local support points
  localIndex numLocalSupport = 0;
  forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
  {
    numLocalSupport += countMeshObjects< false >( field.location, mesh, regions );
  } );
  field.numLocalDof = field.numComponents * numLocalSupport;

  // gather dof counts across ranks
  field.rankOffset = MpiWrapper::prefixSum< globalIndex >( field.numLocalDof );
  field.numGlobalDof = field.rankOffset + field.numLocalDof;
  MpiWrapper::broadcast( field.numGlobalDof, MpiWrapper::commSize() - 1 );

  // determine field's offsets
  if( fieldIndex > 0 )
  {
    FieldDescription const & prev = m_fields[fieldIndex - 1];
    field.blockOffset = prev.blockOffset + prev.numGlobalDof;
  }
  field.globalOffset = field.rankOffset; // actual value computed in reorderByRank()
}

namespace
{

std::vector< string >
processFieldRegionList( MeshLevel const & mesh,
                        std::vector< string > inputList )
{
  std::vector< string > regions( std::move( inputList ) );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  if( regions.empty() )
  {
    elemManager.forElementRegions( [&]( ElementRegionBase const & region )
    {
      regions.push_back( region.getName() );
    } );
  }
  else
  {
    for( string const & regionName : regions )
    {
      // Check existence, discard return value
      elemManager.getRegion( regionName );
    }
  }

  // sort region names and remove duplicates
  std::sort( regions.begin(), regions.end() );
  regions.erase( std::unique( regions.begin(), regions.end() ), regions.end() );

  return regions;
}

std::vector< DofManager::Regions >
processFieldRegionList( DomainPartition const & domain,
                        std::vector< DofManager::Regions > const & inputList )
{
  std::vector< DofManager::Regions > result;
  std::set< std::pair< string, string > > processedMeshLevels;
  for( DofManager::Regions const & r : inputList )
  {
    GEOSX_ERROR_IF( processedMeshLevels.count( { r.meshBodyName, r.meshLevelName } ) > 0,
                    GEOSX_FMT( "Duplicate mesh: {}/{}", r.meshBodyName, r.meshLevelName ) );
    processedMeshLevels.insert( { r.meshBodyName, r.meshLevelName } );

    MeshBody const & body = domain.getMeshBody( r.meshBodyName );
    MeshLevel const & mesh = body.getMeshLevel( r.meshLevelName );
    result.push_back( { r.meshBodyName, r.meshLevelName, processFieldRegionList( mesh, r.regionNames ) } );
  }
  return result;
}

} // namespace

void DofManager::addField( string const & fieldName,
                           Location const location,
                           integer const components,
                           std::vector< Regions > const & regions )
{
  GEOSX_ASSERT_MSG( m_domain != nullptr, "Domain has not been set" );
  GEOSX_ERROR_IF( m_reordered, "Cannot add fields after reorderByRank() has been called." );
  GEOSX_ERROR_IF( fieldExists( fieldName ), "Requested field name '" << fieldName << "' already exists." );
  GEOSX_ERROR_IF_GT_MSG( components, MAX_COMP, "Number of components limit exceeded" );

  m_fields.emplace_back();
  FieldDescription & field = m_fields.back();

  // populate basic info
  field.name = fieldName;
  field.location = location;
  field.numComponents = components;
  field.key = m_name + '_' + fieldName + "_dofIndex";
  field.docstring = fieldName + " DoF indices";

  // advanced processing
  field.support = processFieldRegionList( *m_domain, regions );
  computeFieldDimensions( static_cast< localIndex >( m_fields.size() ) - 1 );

  // allocate and fill index array
  createIndexArray( field );
}

void DofManager::addField( string const & fieldName,
                           Location const location,
                           integer const components,
                           map< string, array1d< string > > const & regions )
{
  // Convert input into internal format
  std::vector< Regions > support;
  for( auto const & p : regions )
  {
    MeshBody const & meshBody = m_domain->getMeshBody( p.first );
    MeshLevel const & mesh = meshBody.getMeshLevel( 0 ); // TODO by name?
    std::vector< string > regionNames( p.second.begin(), p.second.end() );
    support.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
  }
  addField( fieldName, location, components, support );
}

namespace
{

std::vector< string >
processCouplingRegionList( std::vector< string > inputList,
                           std::vector< string > const & rowFieldRegions,
                           string const & rowFieldName,
                           std::vector< string > const & colFieldRegions,
                           string const & colFieldName )
{
  std::vector< string > regions( std::move( inputList ) );
  if( regions.empty() )
  {
    // Populate with regions common between two fields
    regions.resize( std::min( rowFieldRegions.size(), colFieldRegions.size() ) );
    auto const it = std::set_intersection( rowFieldRegions.begin(), rowFieldRegions.end(),
                                           colFieldRegions.begin(), colFieldRegions.end(),
                                           regions.begin() );
    regions.resize( std::distance( regions.begin(), it ) );
  }
  else
  {
    // Sort alphabetically and remove possible duplicates in user input
    std::sort( regions.begin(), regions.end() );
    regions.erase( std::unique( regions.begin(), regions.end() ), regions.end() );

    // Check that both fields exist on all regions in the list
    auto const checkSupport = [&regions]( std::vector< string > const & fieldRegions, string const & fieldName )
    {
      // Both regions lists are sorted at this point
      GEOSX_ERROR_IF( !std::includes( fieldRegions.begin(), fieldRegions.end(), regions.begin(), regions.end() ),
                      GEOSX_FMT( "Coupling domain is not a subset of {}'s support:\nCoupling: {}\nField: {}",
                                 fieldName, stringutilities::join( regions ), stringutilities::join( fieldRegions ) ) );
    };
    checkSupport( rowFieldRegions, rowFieldName );
    checkSupport( colFieldRegions, colFieldName );
  }

  return regions;
}

/// Comparison function that ignores region lists and only looks at mesh body/level names
template< typename OP >
struct RegionComp
{
  bool operator()( DofManager::Regions const & l, DofManager::Regions const & r ) const
  {
    return OP{} ( std::tie( l.meshBodyName, l.meshLevelName ), std::tie( r.meshBodyName, r.meshLevelName ) );
  }
};

std::vector< DofManager::Regions >
processCouplingRegionList( std::vector< DofManager::Regions > inputList,
                           std::vector< DofManager::Regions > const & rowFieldRegions,
                           string const & rowFieldName,
                           std::vector< DofManager::Regions > const & colFieldRegions,
                           string const & colFieldName )
{
  std::vector< DofManager::Regions > regions( std::move( inputList ) );

  if( regions.empty() )
  {
    // First, compute a common set of mesh body/level pairs
    std::set_intersection( rowFieldRegions.begin(), rowFieldRegions.end(),
                           colFieldRegions.begin(), colFieldRegions.end(),
                           std::back_inserter( regions ), RegionComp< std::less<> >{} );
    // Now, compute intersections of region lists
    for( DofManager::Regions & r : regions )
    {
      // Find the body/level pair in the column field list (unsorted range)
      auto const comp = [&r]( auto const & c ){ return RegionComp< std::equal_to<> >{} ( r, c ); };
      DofManager::Regions const & colRegions = *std::find_if( colFieldRegions.begin(), colFieldRegions.end(), comp );
      // Intersect row field regions (already copied into result) with col field regions (found above)
      r.regionNames = processCouplingRegionList( {}, r.regionNames, rowFieldName, colRegions.regionNames, colFieldName );
    }
  }
  else
  {
    // Check that each input entry is included in both row and col field supports
    auto const checkSupport = [&regions]( std::vector< DofManager::Regions > const & fieldRegions, string const & fieldName )
    {
      for( DofManager::Regions const & r : regions )
      {
        auto const comp = [&r]( auto const & f ){ return RegionComp< std::equal_to<> >{} ( r, f ); };
        auto const it = std::find_if( fieldRegions.begin(), fieldRegions.end(), comp );
        GEOSX_ERROR_IF( it == fieldRegions.end(),
                        GEOSX_FMT( "Mesh {}/{} not found in support of field {}", r.meshBodyName, r.meshLevelName, fieldName ) );
        GEOSX_ERROR_IF( !std::includes( it->regionNames.begin(), it->regionNames.end(), r.regionNames.begin(), r.regionNames.end() ),
                        GEOSX_FMT( "Coupling domain is not a subset of {}'s support:\nCoupling: {}\nField: {}",
                                   fieldName, stringutilities::join( r.regionNames ), stringutilities::join( it->regionNames ) ) );
      }
    };
    checkSupport( rowFieldRegions, rowFieldName );
    checkSupport( colFieldRegions, colFieldName );
  }

  return regions;
}

} // namespace

void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity,
                              std::vector< Regions > const & regions,
                              bool const symmetric )
{
  GEOSX_ASSERT_MSG( m_domain != nullptr, "Mesh has not been set" );
  localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
  localIndex const colFieldIndex = getFieldIndex( colFieldName );

  // Check if already defined
  GEOSX_ASSERT_MSG( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0, "addCoupling: coupling already defined" );

  CouplingDescription & coupling = m_coupling[ { rowFieldIndex, colFieldIndex } ];
  coupling.connector = connectivity;
  if( connectivity == Connector::None && rowFieldIndex == colFieldIndex )
  {
    coupling.connector = static_cast< Connector >( m_fields[rowFieldIndex].location );
  }

  // Make a list of regions on which coupling is defined
  coupling.support = processCouplingRegionList( regions,
                                                m_fields[rowFieldIndex].support,
                                                rowFieldName,
                                                m_fields[colFieldIndex].support,
                                                colFieldName );

  // Set connectivity with active symmetry flag
  if( symmetric && colFieldIndex != rowFieldIndex )
  {
    m_coupling.insert( { { colFieldIndex, rowFieldIndex }, coupling } );
  }
}

void DofManager::addCoupling( string const & fieldName,
                              FluxApproximationBase const & stencils )
{
  localIndex const fieldIndex = getFieldIndex( fieldName );
  FieldDescription const & field = m_fields[fieldIndex];

  GEOSX_ERROR_IF( field.location != Location::Elem, "Field must be supported on elements in order to use stencil sparsity" );

  CouplingDescription & coupling = m_coupling[ {fieldIndex, fieldIndex} ];
  coupling.connector = Connector::Stencil;
  coupling.support = field.support;
  coupling.stencils = &stencils;
}

void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              DofManager::Connector connectivity,
                              map< string, array1d< string > > const & regions,
                              bool symmetric )
{
  // Convert input into internal format
  std::vector< Regions > support;
  for( auto const & p : regions )
  {
    MeshBody const & meshBody = m_domain->getMeshBody( p.first );
    MeshLevel const & mesh = meshBody.getMeshLevel( 0 ); // TODO by name?
    std::vector< string > regionNames( p.second.begin(), p.second.end() );
    support.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
  }
  addCoupling( rowFieldName, colFieldName, connectivity, support, symmetric );
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
  static void build( MeshLevel const & mesh,
                     string const & key,
                     localIndex const numComp,
                     std::vector< string > const & regions,
                     localIndex const rowOffset,
                     SparsityPattern< globalIndex > & connLocPattern )
  {
    using helper = ArrayHelper< globalIndex const, LOC >;
    typename helper::Accessor dofIndexArray = helper::get( mesh, key );
    auto connIdxPrev = MeshHelper< CONN >::invalid_local_index;

    localIndex connectorCount = -1;

    forMeshLocation< CONN, LOC, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions,
                                                                         [&]( auto const & connIdx,
                                                                              auto const & locIdx,
                                                                              localIndex const GEOSX_UNUSED_PARAM( k ) )
    {
      if( connIdx != connIdxPrev )
      {
        ++connectorCount;
        connIdxPrev = connIdx;
      }

      globalIndex const dofOffset = helper::value( dofIndexArray, locIdx );
      if( dofOffset >= 0 )
      {
        for( localIndex c = 0; c < numComp; ++c )
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
  static void build( MeshLevel const & mesh,
                     string const & key,
                     localIndex const numComp,
                     std::vector< string > const & regions,
                     localIndex const rowOffset,
                     SparsityPattern< globalIndex > & connLocPattern )
  {
    DofManager::Location constexpr ELEM  = DofManager::Location::Elem;
    DofManager::Location constexpr EDGE = DofManager::Location::Edge;

    EdgeManager const & edgeManager = mesh.getEdgeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofIndex =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( key );

    array1d< localIndex > edgeConnectorIndex( edgeManager.size() );
    edgeConnectorIndex.setValues< serialPolicy >( -1 );

    localIndex edgeCount = 0;
    forMeshLocation< EDGE, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions, [&] ( localIndex const edgeIdx )
    {
      edgeConnectorIndex[edgeIdx] = edgeCount++;
    } );

    forMeshLocation< ELEM, EDGE, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions,
                                                                          [&]( auto const & elemIdx,
                                                                               auto const & edgeIdx,
                                                                               localIndex const GEOSX_UNUSED_PARAM( k ) )
    {
      globalIndex const dofOffset = dofIndex[elemIdx[0]][elemIdx[1]][elemIdx[2]];
      if( dofOffset >= 0 )
      {
        for( localIndex c = 0; c < numComp; ++c )
        {
          connLocPattern.insertNonZero( rowOffset + edgeConnectorIndex[edgeIdx], dofOffset + c );
        }
      }
    } );
  }
};

template< DofManager::Location LOC, DofManager::Location CONN >
void makeConnLocPattern( MeshLevel const & mesh,
                         string const & key,
                         localIndex const numComp,
                         globalIndex const numGlobalDof,
                         std::vector< string > const & regions,
                         SparsityPattern< globalIndex > & connLocPattern )
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
  localIndex const numEntriesPerRow = MeshIncidence< CONN, LOC >::max * numComp;
  connLocPattern.resize( numConnectorsTotal,
                         numGlobalDof,
                         std::min( LvArray::integerConversion< globalIndex >( numEntriesPerRow ), numGlobalDof ) );

  // 3. Populate the local CL pattern
  ConnLocPatternBuilder< LOC, CONN >::build( mesh, key, numComp, regions, 0, connLocPattern );

  // TPFA+fracture special case
  if( LOC == Loc::Elem && CONN == Loc::Face )
  {
    ConnLocPatternBuilder< Loc::Elem, Loc::Edge, FaceElementSubRegion >::build( mesh, key, numComp, regions, numConnectors, connLocPattern );
  }
}

} // namespace

void DofManager::setSparsityPatternFromStencil( SparsityPatternView< globalIndex > const & pattern,
                                                localIndex const fieldIndex ) const
{
  FieldDescription const & field = m_fields[fieldIndex];
  CouplingDescription const & coupling = m_coupling.at( {fieldIndex, fieldIndex} );
  localIndex const numComp = field.numComponents;
  globalIndex const rankDofOffset = rankOffset();

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const dofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( field.key );

    array1d< globalIndex > rowDofIndices( numComp );
    array1d< globalIndex > colDofIndices( numComp );

    // 1. Assemble diagonal and off-diagonal blocks for elements in stencil
    coupling.stencils->forAllStencils( mesh, [&]( auto const & stencil )
    {
      using StenciType = typename std::decay< decltype( stencil ) >::type;
      constexpr localIndex maxNumFluxElems = StenciType::NUM_POINT_IN_FLUX;
      constexpr localIndex maxStencilSize = StenciType::MAX_STENCIL_SIZE;

      typename StenciType::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      rowDofIndices.reserve( maxNumFluxElems );
      colDofIndices.reserve( maxStencilSize * numComp );

      forAll< serialPolicy >( stencil.size(), [&]( localIndex const iconn )
      {
        // This weirdness is because of fracture stencils, which don't have separate
        // getters for num flux elems vs stencil size... it won't work for MPFA though
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        localIndex const stencilSize = numFluxElems;

        rowDofIndices.resize( numFluxElems );
        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          rowDofIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        }

        colDofIndices.resize( stencilSize * numComp );
        for( localIndex i = 0; i < stencilSize; ++i )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            colDofIndices[i * numComp + c] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + c;
          }
        }

        std::sort( colDofIndices.begin(), colDofIndices.end() );

        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          localIndex const localDofNumber = rowDofIndices[i] - rankDofOffset;
          if( localDofNumber >= 0 && localDofNumber < pattern.numRows() )
          {
            for( localIndex c = 0; c < numComp; ++c )
            {
              pattern.insertNonZeros( localDofNumber + c, colDofIndices.begin(), colDofIndices.end() );
            }
          }
        }
      } );
    } );

    // 2. Insert diagonal blocks, in case there are elements not included in stencil
    // (e.g. a single fracture element not connected to any other)
    auto dofNumberView = dofNumber.toNestedViewConst();
    colDofIndices.resize( numComp );
    forMeshLocation< Location::Elem, false, parallelHostPolicy >( mesh, regions, [=]( auto const & elemIdx )
    {
      globalIndex const elemDof = dofNumberView[elemIdx[0]][elemIdx[1]][elemIdx[2]];
      for( localIndex c = 0; c < numComp; ++c )
      {
        colDofIndices[c] = elemDof + c;
      }
      for( localIndex c = 0; c < numComp; ++c )
      {
        pattern.insertNonZeros( elemDof - rankDofOffset + c, colDofIndices.begin(), colDofIndices.end() );
      }
    } );
  } );
}

void DofManager::setSparsityPatternOneBlock( SparsityPatternView< globalIndex > const & pattern,
                                             localIndex const rowFieldIndex,
                                             localIndex const colFieldIndex ) const
{
  GEOSX_ASSERT( rowFieldIndex >= 0 );
  GEOSX_ASSERT( colFieldIndex >= 0 );

  FieldDescription const & rowField = m_fields[rowFieldIndex];
  FieldDescription const & colField = m_fields[colFieldIndex];

  if( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0 )
  {
    return;
  }
  CouplingDescription const & coupling = m_coupling.at( {rowFieldIndex, colFieldIndex} );

  // Special treatment for stencil-based sparsity
  if( rowFieldIndex == colFieldIndex && coupling.connector == Connector::Stencil )
  {
    setSparsityPatternFromStencil( pattern, rowFieldIndex );
    return;
  }

  forMeshSupport( coupling.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    SparsityPattern< globalIndex > connLocRow( 0, 0, 0 ), connLocCol( 0, 0, 0 );

    LocationSwitch( rowField.location, static_cast< Location >( coupling.connector ),
                    [&]( auto const locType, auto const connType )
    {
      Location constexpr LOC  = decltype(locType)::value;
      Location constexpr CONN = decltype(connType)::value;

      makeConnLocPattern< LOC, CONN >( mesh,
                                       rowField.key,
                                       rowField.numComponents,
                                       rowField.numGlobalDof,
                                       regions,
                                       connLocRow );
    } );

    if( colFieldIndex == rowFieldIndex )
    {
      connLocCol = connLocRow; // TODO avoid copying
    }
    else
    {
      LocationSwitch( colField.location, static_cast< Location >( coupling.connector ),
                      [&]( auto const locType, auto const connType )
      {
        Location constexpr LOC = decltype(locType)::value;
        Location constexpr CONN = decltype(connType)::value;

        makeConnLocPattern< LOC, CONN >( mesh,
                                         colField.key,
                                         colField.numComponents,
                                         colField.numGlobalDof,
                                         regions,
                                         connLocCol );
      } );
    }
    GEOSX_ASSERT_EQ( connLocRow.numRows(), connLocCol.numRows() );

    // Perform assembly/multiply patterns
    globalIndex const globalDofOffset = rankOffset();
    forAll< serialPolicy >( connLocRow.numRows(), [&]( localIndex const irow )
    {
      arraySlice1d< globalIndex const > const dofIndicesRow = connLocRow.getColumns( irow );
      arraySlice1d< globalIndex const > const dofIndicesCol = connLocCol.getColumns( irow );
      for( globalIndex const globalRow : dofIndicesRow )
      {
        localIndex const localRow = globalRow - globalDofOffset;
        if( localRow >= 0 && localRow < pattern.numRows() )
        {
          pattern.insertNonZeros( localRow, dofIndicesCol.begin(), dofIndicesCol.end() );
        }
      }
    } );
  } );
}

void DofManager::countRowLengthsFromStencil( arrayView1d< localIndex > const & rowLengths,
                                             localIndex const fieldIndex ) const
{
  FieldDescription const & field = m_fields[fieldIndex];
  CouplingDescription const & coupling = m_coupling.at( {fieldIndex, fieldIndex} );
  GEOSX_ASSERT( coupling.connector == Connector::Stencil );
  GEOSX_ASSERT( coupling.stencils != nullptr );

  localIndex const numComp = field.numComponents;
  globalIndex const rankDofOffset = rankOffset();

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const dofNumberAccessor =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( field.key );

    // 1. Count row contributions from stencil
    coupling.stencils->forAllStencils( mesh, [&]( auto const & stencil )
    {
      using StenciType = typename std::decay< decltype( stencil ) >::type;
      typename StenciType::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      forAll< parallelHostPolicy >( stencil.size(), [&]( localIndex const iconn )
      {
        // This weirdness is because of fracture stencils, which don't have separate
        // getters for num flux elems vs stencil size... it won't work for MPFA though
        localIndex const stencilSize = stencil.stencilSize( iconn );
        localIndex const numFluxElems = stencilSize;

        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          localIndex const localDofNumber = dofNumberAccessor[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] - rankDofOffset;
          if( localDofNumber >= 0 && localDofNumber < rowLengths.size() )
          {
            for( localIndex c = 0; c < numComp; ++c )
            {
              RAJA::atomicAdd( parallelHostAtomic{}, &rowLengths[localDofNumber + c], ( stencilSize - 1 ) * numComp );
            }
          }
        }
      } );
    } );

    // 2. Add diagonal contributions to account for elements not in stencil
    mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const, ElementSubRegionBase const & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( field.key );
      forAll< parallelHostPolicy >( subRegion.size(), [&]( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          localIndex const localDofNumber = dofNumber[ei] - rankDofOffset;
          for( localIndex c = 0; c < numComp; ++c )
          {
            rowLengths[localDofNumber + c] += numComp;
          }
        }
      } );
    } );
  } );
}

void DofManager::countRowLengthsOneBlock( arrayView1d< localIndex > const & rowLengths,
                                          localIndex const rowFieldIndex,
                                          localIndex const colFieldIndex ) const
{
  GEOSX_ASSERT( rowFieldIndex >= 0 );
  GEOSX_ASSERT( colFieldIndex >= 0 );

  FieldDescription const & rowField = m_fields[rowFieldIndex];
  FieldDescription const & colField = m_fields[colFieldIndex];

  if( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0 )
  {
    return;
  }
  CouplingDescription const & coupling = m_coupling.at( {rowFieldIndex, colFieldIndex} );

  if( rowFieldIndex == colFieldIndex && coupling.connector == Connector::Stencil )
  {
    countRowLengthsFromStencil( rowLengths, rowFieldIndex );
    return;
  }

  forMeshSupport( coupling.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    SparsityPattern< globalIndex > connLocRow( 0, 0, 0 ), connLocCol( 0, 0, 0 );

    LocationSwitch( rowField.location, static_cast< Location >( coupling.connector ),
                    [&]( auto const locType, auto const connType )
    {
      Location constexpr LOC  = decltype(locType)::value;
      Location constexpr CONN = decltype(connType)::value;

      makeConnLocPattern< LOC, CONN >( mesh,
                                       rowField.key,
                                       rowField.numComponents,
                                       rowField.numGlobalDof,
                                       regions,
                                       connLocRow );
    } );

    if( colFieldIndex == rowFieldIndex )
    {
      connLocCol = connLocRow;
    }
    else
    {
      LocationSwitch( colField.location, static_cast< Location >( coupling.connector ),
                      [&]( auto const locType, auto const connType )
      {
        Location constexpr LOC = decltype(locType)::value;
        Location constexpr CONN = decltype(connType)::value;

        makeConnLocPattern< LOC, CONN >( mesh,
                                         colField.key,
                                         colField.numComponents,
                                         colField.numGlobalDof,
                                         regions,
                                         connLocCol );
      } );
    }

    GEOSX_ASSERT_EQ( connLocRow.numRows(), connLocCol.numRows() );

    // Estimate an upper bound on row length by adding adjacent entries on a connector
    globalIndex const rankDofOffset = rankOffset();
    forAll< parallelHostPolicy >( connLocRow.numRows(), [&]( localIndex const iconn )
    {
      arraySlice1d< globalIndex const > const dofIndicesRow = connLocRow.getColumns( iconn );
      for( globalIndex const globalRow : dofIndicesRow )
      {
        localIndex const localRow = globalRow - rankDofOffset;
        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          RAJA::atomicAdd( parallelHostAtomic{}, &rowLengths[localRow], connLocCol.numNonZeros( iconn ) );
        }
      }
    } );
  } );
}

// Create the sparsity pattern (location-location). Low level interface
void DofManager::setSparsityPattern( SparsityPattern< globalIndex > & pattern ) const
{
  GEOSX_ERROR_IF( !m_reordered, "Cannot set monolithic sparsity pattern before reorderByRank() has been called." );

  localIndex const numLocalRows = numLocalDofs();
  localIndex const numFields = LvArray::integerConversion< localIndex >( m_fields.size() );

  // Step 1. Do a dry run of sparsity construction to get the total number of nonzeros in each row
  array1d< localIndex > rowSizes( numLocalRows );
  for( localIndex blockRow = 0; blockRow < numFields; ++blockRow )
  {
    for( localIndex blockCol = 0; blockCol < numFields; ++blockCol )
    {
      countRowLengthsOneBlock( rowSizes, blockRow, blockCol );
    }
  }

  // Step 2. Allocate enough capacity for all nonzero entries in each row
  pattern.resizeFromRowCapacities< parallelHostPolicy >( numLocalRows, numGlobalDofs(), rowSizes.data() );

  // Step 3. Fill the sparsity block-by-block
  for( localIndex blockRow = 0; blockRow < numFields; ++blockRow )
  {
    for( localIndex blockCol = 0; blockCol < numFields; ++blockCol )
    {
      setSparsityPatternOneBlock( pattern.toView(), blockRow, blockCol );
    }
  }

  // Step 4. Compress to remove unused space between rows
  pattern.compress();
}

namespace
{

template< typename FIELD_OP, typename POLICY, typename FIELD_VIEW >
void vectorToFieldKernel( arrayView1d< real64 const > const & localVector,
                          FIELD_VIEW const & field,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< integer const > const & ghostRank,
                          real64 const scalingFactor,
                          localIndex const dofOffset,
                          DofManager::CompMask const mask )
{
  forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    if( ghostRank[i] < 0 && dofNumber[i] >= 0 )
    {
      localIndex const lid = dofNumber[i] - dofOffset;
      GEOSX_ASSERT_GE( lid, 0 );

      integer fieldComp = 0;
      for( integer const vecComp : mask )
      {
        FIELD_OP::template SpecifyFieldValue( field,
                                              i,
                                              fieldComp++,
                                              scalingFactor * localVector[lid + vecComp] );
      }
    }
  } );
}

template< typename FIELD_OP, typename POLICY >
void vectorToFieldImpl( arrayView1d< real64 const > const & localVector,
                        ObjectManagerBase & manager,
                        string const & dofKey,
                        string const & fieldName,
                        real64 const scalingFactor,
                        localIndex const dofOffset,
                        DofManager::CompMask const mask )
{
  arrayView1d< globalIndex const > const dofNumber = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const ghostRank = manager.ghostRank();

  WrapperBase & wrapper = manager.getWrapperBase( fieldName );

  // Restrict primary solution fields to 1-2D real arrays,
  // because applying component index is not well defined for 3D and higher
  using FieldTypes = types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > >;
  types::dispatch( FieldTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    vectorToFieldKernel< FIELD_OP, POLICY >( localVector,
                                             wrapperT.reference().toView(),
                                             dofNumber,
                                             ghostRank,
                                             scalingFactor,
                                             dofOffset,
                                             mask );
  } );
}

template< typename FIELD_OP, typename POLICY, typename FIELD_VIEW >
void fieldToVectorKernel( arrayView1d< real64 > const & localVector,
                          FIELD_VIEW const & field,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< integer const > const & ghostRank,
                          real64 const GEOSX_UNUSED_PARAM( scalingFactor ),
                          localIndex const dofOffset,
                          DofManager::CompMask const mask )
{
  forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    if( ghostRank[i] < 0 && dofNumber[i] >= 0 )
    {
      localIndex const lid = dofNumber[i] - dofOffset;
      GEOSX_ASSERT_GE( lid, 0 );

      integer fieldComp = 0;
      for( integer const vecComp : mask )
      {
        FIELD_OP::template readFieldValue( field,
                                           i,
                                           fieldComp++,
                                           localVector[lid + vecComp] );
      }
    }
  } );
}

template< typename FIELD_OP, typename POLICY >
void fieldToVectorImpl( arrayView1d< real64 > const & localVector,
                        ObjectManagerBase const & manager,
                        string const & dofKey,
                        string const & fieldName,
                        real64 const scalingFactor,
                        localIndex const dofOffset,
                        DofManager::CompMask const mask )
{
  arrayView1d< globalIndex const > const & dofNumber = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = manager.ghostRank();

  WrapperBase const & wrapper = manager.getWrapperBase( fieldName );

  // Restrict primary solution fields to 1-2D real arrays,
  // because applying component index is not well defined for 3D and higher
  using FieldTypes = types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > >;
  types::dispatch( FieldTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > const & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    fieldToVectorKernel< FIELD_OP, POLICY >( localVector,
                                             wrapperT.reference(),
                                             dofNumber,
                                             ghostRank,
                                             scalingFactor,
                                             dofOffset,
                                             mask );
  } );
}

} // namespace

template< typename FIELD_OP, typename POLICY >
void DofManager::vectorToField( arrayView1d< real64 const > const & localVector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                real64 const scalingFactor,
                                CompMask mask ) const
{
  FieldDescription const & field = m_fields[getFieldIndex( srcFieldName )];
  mask.setNumComp( field.numComponents );

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel & mesh, auto const & regions )
  {
    if( field.location == Location::Elem )
    {
      mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
      {
        vectorToFieldImpl< FIELD_OP, POLICY >( localVector,
                                               subRegion,
                                               field.key,
                                               dstFieldName,
                                               scalingFactor,
                                               rankOffset(),
                                               mask );
      } );
    }
    else
    {
      vectorToFieldImpl< FIELD_OP, POLICY >( localVector,
                                             getObjectManager( field.location, mesh ),
                                             field.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             mask );
    }
  } );
}

void DofManager::copyVectorToField( arrayView1d< real64 const > const & localVector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    CompMask const mask ) const
{
  vectorToField< FieldSpecificationEqual, parallelDevicePolicy<> >( localVector,
                                                                    srcFieldName,
                                                                    dstFieldName,
                                                                    scalingFactor,
                                                                    mask );
}

void DofManager::addVectorToField( arrayView1d< real64 const > const & localVector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   real64 const scalingFactor,
                                   CompMask const mask ) const
{
  vectorToField< FieldSpecificationAdd, parallelDevicePolicy<> >( localVector,
                                                                  srcFieldName,
                                                                  dstFieldName,
                                                                  scalingFactor,
                                                                  mask );
}

template< typename FIELD_OP, typename POLICY >
void DofManager::fieldToVector( arrayView1d< real64 > const & localVector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                real64 const scalingFactor,
                                CompMask mask ) const
{
  FieldDescription const & field = m_fields[getFieldIndex( srcFieldName )];
  mask.setNumComp( field.numComponents );

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel & mesh, auto const & regions )
  {
    if( field.location == Location::Elem )
    {
      mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const,
                                                                ElementSubRegionBase const & subRegion )
      {
        fieldToVectorImpl< FIELD_OP, POLICY >( localVector,
                                               subRegion,
                                               field.key,
                                               dstFieldName,
                                               scalingFactor,
                                               rankOffset(),
                                               mask );
      } );
    }
    else
    {
      fieldToVectorImpl< FIELD_OP, POLICY >( localVector,
                                             getObjectManager( field.location, mesh ),
                                             field.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             mask );
    }
  } );
}

void DofManager::copyFieldToVector( arrayView1d< real64 > const & localVector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    CompMask const mask ) const
{
  fieldToVector< FieldSpecificationEqual, parallelDevicePolicy<> >( localVector,
                                                                    srcFieldName,
                                                                    dstFieldName,
                                                                    scalingFactor,
                                                                    mask );
}

void DofManager::addFieldToVector( arrayView1d< real64 > const & localVector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   real64 const scalingFactor,
                                   CompMask const mask ) const
{
  fieldToVector< FieldSpecificationAdd, parallelDevicePolicy<> >( localVector,
                                                                  srcFieldName,
                                                                  dstFieldName,
                                                                  scalingFactor,
                                                                  mask );
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

  std::map< std::pair< string, string >, std::map< string, string_array > > fieldsToSync;

  // adjust index arrays for owned locations
  for( FieldDescription const & field : m_fields )
  {
    globalIndex const adjustment = field.globalOffset - field.rankOffset;

    forMeshSupport( field.support, *m_domain, [&]( MeshBody const & body, MeshLevel & mesh, auto const & regions )
    {
      LocationSwitch( field.location, [&]( auto const loc )
      {
        Location constexpr LOC = decltype(loc)::value;
        using ArrayHelper = ArrayHelper< globalIndex, LOC >;

        typename ArrayHelper::Accessor indexArray = ArrayHelper::get( mesh, field.key );

        forMeshLocation< LOC, false, parallelHostPolicy >( mesh, regions, [&]( auto const locIdx )
        {
          ArrayHelper::reference( indexArray, locIdx ) += adjustment;
        } );

        fieldsToSync[{ body.getName(), mesh.getName() }][MeshHelper< LOC >::syncObjName].emplace_back( field.key );
      } );
    } );
  }

  // synchronize index arrays for all fields across ranks
  for( auto const & meshFieldPair : fieldsToSync )
  {
    MeshLevel & mesh = m_domain->getMeshBody( meshFieldPair.first.first ).getMeshLevel( meshFieldPair.first.second );
    CommunicationTools::getInstance().synchronizeFields( meshFieldPair.second, mesh, m_domain->getNeighbors(), false );
  }

  m_reordered = true;
}

std::vector< DofManager::SubComponent >
DofManager::filterDofs( std::vector< SubComponent > const & excluded ) const
{
  std::vector< DofManager::SubComponent > result;
  for( const auto & field : m_fields )
  {
    auto const it = std::find_if( excluded.begin(), excluded.end(),
                                  [&]( SubComponent const & sc ) { return sc.fieldName == field.name; } );
    CompMask const includedComp = ( it != excluded.end() ) ? ~it->mask : ~CompMask( field.numComponents );
    if( !includedComp.empty() )
    {
      result.push_back( { field.name, includedComp } );
    }
  }

  return result;
}

void DofManager::setupFrom( DofManager const & source,
                            std::vector< SubComponent > const & selection )
{
  clear();
  for( FieldDescription const & field : source.m_fields )
  {
    auto const it = std::find_if( selection.begin(), selection.end(),
                                  [&]( SubComponent const & sc ) { return sc.fieldName == field.name; } );
    if( it != selection.end() )
    {
      SubComponent const & sc = *it;
      GEOSX_ASSERT_GE( field.numComponents, sc.mask.size() );
      addField( field.name, field.location, sc.mask.size(), field.support );
    }
  }
  for( auto const & entry : source.m_coupling )
  {
    string const & rowFieldName = source.m_fields[entry.first.first].name;
    string const & colFieldName = source.m_fields[entry.first.second].name;
    if( fieldExists( rowFieldName ) && fieldExists( colFieldName ) )
    {
      localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
      localIndex const colFieldIndex = getFieldIndex( colFieldName );
      m_coupling.insert( { {rowFieldIndex, colFieldIndex}, entry.second } );
    }
  }
  reorderByRank();
}

template< typename MATRIX >
void DofManager::makeRestrictor( std::vector< SubComponent > const & selection,
                                 MPI_Comm const & comm,
                                 bool const transpose,
                                 MATRIX & restrictor ) const
{
  GEOSX_ERROR_IF( !m_reordered, "Cannot make restrictors before reorderByRank() has been called." );

  // 1. Populate selected fields and compute some basic dimensions
  // array1d< FieldDescription > fieldsSelected( selection.size() );
  std::vector< FieldDescription > fieldsSelected( selection.size() );

  for( std::size_t k = 0; k < fieldsSelected.size(); ++k )
  {
    SubComponent const & dof = selection[k];
    FieldDescription const & fieldOld = m_fields[getFieldIndex( dof.fieldName )];
    FieldDescription & fieldNew = fieldsSelected[k];

    fieldNew.name = selection[k].fieldName;
    fieldNew.numComponents = dof.mask.size();
    fieldNew.numLocalDof = fieldOld.numLocalDof / fieldOld.numComponents * fieldNew.numComponents;
    fieldNew.numGlobalDof = fieldOld.numGlobalDof / fieldOld.numComponents * fieldNew.numComponents;
    fieldNew.rankOffset = fieldOld.rankOffset / fieldOld.numComponents * fieldNew.numComponents;
  }

  // 2. Compute remaining offsets (we mostly just need globalOffset, but it depends on others)

  globalIndex blockOffset = 0;
  for( auto & field : fieldsSelected )
  {
    field.blockOffset = blockOffset;
    blockOffset += field.numGlobalDof;
  }

  globalIndex globalOffset = std::accumulate( fieldsSelected.begin(), fieldsSelected.end(), globalIndex( 0 ),
                                              []( localIndex const n, FieldDescription const & f )
  { return n + f.rankOffset; } );

  for( auto & field : fieldsSelected )
  {
    field.globalOffset = globalOffset;
    globalOffset += field.numLocalDof;
  }

  // 3. Build the restrictor field by field

  localIndex const numLocalDofSelected = std::accumulate( fieldsSelected.begin(), fieldsSelected.end(), localIndex( 0 ),
                                                          []( localIndex const n, FieldDescription const & f )
  { return n + f.numLocalDof; } );

  localIndex const rowSize = transpose ? numLocalDofs() : numLocalDofSelected;
  localIndex const colSize = transpose ? numLocalDofSelected : numLocalDofs();

  restrictor.createWithLocalSize( rowSize, colSize, 1, comm );
  restrictor.open();

  array1d< globalIndex > rows;
  array1d< globalIndex > cols;
  array1d< real64 > values;

  for( std::size_t k = 0; k < fieldsSelected.size(); ++k )
  {
    FieldDescription const & fieldNew = fieldsSelected[k];
    FieldDescription const & fieldOld = m_fields[getFieldIndex( fieldNew.name )];

    FieldDescription const & fieldRow = transpose ? fieldOld : fieldNew;
    FieldDescription const & fieldCol = transpose ? fieldNew : fieldOld;

    CompMask const mask = selection[k].mask;
    localIndex const numLocalNodes = fieldNew.numLocalDof / fieldNew.numComponents;


    rows.resize( numLocalNodes*mask.size() );
    cols.resize( numLocalNodes*mask.size() );
    values.resize( numLocalNodes*mask.size() );

    for( localIndex i = 0; i < numLocalNodes; ++i )
    {
      integer newComp = 0;
      for( integer const oldComp : mask )
      {
        globalIndex const row = fieldRow.globalOffset + i * fieldRow.numComponents + ( transpose ? oldComp : newComp );
        globalIndex const col = fieldCol.globalOffset + i * fieldCol.numComponents + ( transpose ? newComp : oldComp );

        rows( i*mask.size() + newComp ) = row;
        cols( i*mask.size() + newComp ) = col;
        values[ i*mask.size() + newComp ] = 1.0;
        ++newComp;
      }
    }

    restrictor.insert( rows.toViewConst(),
                       cols.toViewConst(),
                       values.toViewConst() );


  }
  restrictor.close();
}


void DofManager::printFieldInfo( std::ostream & os ) const
{
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    localIndex const numFields = LvArray::integerConversion< localIndex >( m_fields.size() );

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
        if( m_coupling.count( {i, j} ) == 0 )
        {
          continue;
        }
        switch( m_coupling.at( {i, j} ).connector )
        {
          case Connector::Elem:
          {
            os << " E ";
            break;
          }
          case Connector::Face:
          {
            os << " F ";
            break;
          }
          case Connector::Edge:
          {
            os << " G ";
            break;
          }
          case Connector::Node:
          {
            os << " N ";
            break;
          }
          case Connector::None:
          {
            os << "   ";
            break;
          }
          case Connector::Stencil:
          {
            os << " S ";
            break;
          }
          default:
          {
            GEOSX_ERROR( "Invalid connector type: " << static_cast< int >( m_coupling.at( {i, j} ).connector ) );
          }
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
  template void DofManager::makeRestrictor( std::vector< SubComponent > const & selection, \
                                            MPI_Comm const & comm, \
                                            bool const transpose, \
                                            LAI::ParallelMatrix & restrictor ) const;

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
