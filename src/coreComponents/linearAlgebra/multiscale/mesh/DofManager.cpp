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
 * @file DogManager.hpp
 */

#include "DofManager.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "MeshObjectManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{
namespace multiscale
{

DofManager::DofManager( string name )
  : m_name( std::move( name ) )
{}

void DofManager::clear()
{
  // deallocate index arrays from the mesh
  for( FieldDescription const & field : m_fields )
  {
    field.manager->deregisterWrapper( field.key );
  }

  // delete internal data
  m_fields.clear();
  m_reordered = false;
}

void DofManager::setDomain( DomainPartition & domain )
{
  clear();
  m_domain = &domain;
}

DofManager::FieldDescription const & DofManager::getField( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  GEOS_ASSERT_MSG( it != m_fields.end(), "DofManager: field does not exist: " << name );
  return *it;
}

bool DofManager::fieldExists( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  return it != m_fields.end();
}

string const & DofManager::key( string const & fieldName ) const
{
  return getField( fieldName ).key;
}

globalIndex DofManager::numGlobalDofs( string const & fieldName ) const
{
  return getField( fieldName ).numGlobalDof;
}

globalIndex DofManager::numGlobalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex{},
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numGlobalDof; } );
}

localIndex DofManager::numLocalDofs( string const & fieldName ) const
{
  return getField( fieldName ).numLocalDof;
}

localIndex DofManager::numLocalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), localIndex{},
                          []( localIndex const n, FieldDescription const & f ) { return n + f.numLocalDof; } );
}

localIndex DofManager::localOffset( string const & fieldName ) const
{
  return getField( fieldName ).localOffset;
}

globalIndex DofManager::rankOffset( string const & fieldName ) const
{
  return getField( fieldName ).rankOffset;
}

globalIndex DofManager::rankOffset() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex{},
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.rankOffset; } );
}

integer DofManager::numComponents( string const & fieldName ) const
{
  return getField( fieldName ).numComponents;
}

integer DofManager::numComponents() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), integer{},
                          []( integer const n, FieldDescription const & f ) { return n + f.numComponents; } );
}

globalIndex DofManager::globalOffset( string const & fieldName ) const
{
  GEOS_ASSERT_MSG( m_reordered, "Global offset not available until after reorderByRank() has been called." );
  return getField( fieldName ).globalOffset;
}

MeshObjectManager const & DofManager::manager( string const & fieldName ) const
{
  return *getField( fieldName ).manager;
}

void DofManager::createIndexArray( FieldDescription const & field )
{
  // register index array
  arrayView1d< globalIndex > const indexArray =
    field.manager->registerWrapper< array1d< globalIndex > >( field.key ).
      setApplyDefaultValue( -1 ).
      setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
      setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
      setDescription( "Degree-of-freedom indices for " + field.name ).
      referenceAsView();

  // Assign local DoF indices
  forAll< parallelHostPolicy >( field.manager->numOwnedObjects(), [indexArray, &field]( localIndex const i )
  {
    indexArray[i] = field.rankOffset + field.numComponents * i;
  } );

  // synchronize across ranks
  string_array fieldNames;
  fieldNames.emplace_back( field.key );
  CommunicationTools::getInstance().synchronizeFields( fieldNames, *field.manager, m_domain->getNeighbors(), false );
}

void DofManager::addField( string const & fieldName,
                           integer const components,
                           MeshObjectManager & manager )
{
  GEOS_ERROR_IF( m_domain == nullptr, "Domain has not been set" );
  GEOS_ERROR_IF( m_reordered, "Cannot add fields after reorderByRank() has been called." );
  GEOS_ERROR_IF( fieldExists( fieldName ), GEOS_FMT( "Duplicate field declared: '{}'.", fieldName ) );
  GEOS_ERROR_IF_LT_MSG( components, 1, "Invalid number of components specified." );

  m_fields.emplace_back();
  FieldDescription & field = m_fields.back();

  // Fill basic information
  field.name = fieldName;
  field.manager = &manager;
  field.numComponents = components;
  field.key = m_name + '_' + fieldName + "_dofIndex";

  // determine number of local support points
  localIndex const numLocalSupport = field.manager->numOwnedObjects();
  field.numLocalDof = field.numComponents * numLocalSupport;

  // gather dof counts across ranks
  field.rankOffset = MpiWrapper::prefixSum< globalIndex >( field.numLocalDof );
  field.numGlobalDof = field.rankOffset + field.numLocalDof;
  MpiWrapper::broadcast( field.numGlobalDof, MpiWrapper::commSize() - 1 );

  // determine field's offsets
  if( m_fields.size() > 1 )
  {
    FieldDescription const & prev = m_fields[m_fields.size() - 2];
    field.localOffset = prev.localOffset + prev.numLocalDof;
    field.blockOffset = prev.blockOffset + prev.numGlobalDof;
  }
  field.globalOffset = field.rankOffset; // actual value computed in reorderByRank()

  createIndexArray( field );
}

void DofManager::reorderByRank()
{
  GEOS_ASSERT( !m_reordered );

  m_reordered = true;
  if( m_fields.size() < 2 )
  {
    return;
  }

  // update field offsets to account for renumbering
  globalIndex dofOffset = rankOffset();
  for( FieldDescription & field : m_fields )
  {
    field.globalOffset = dofOffset;
    dofOffset += field.numLocalDof;
  }

  // adjust index arrays for owned locations
  for( FieldDescription const & field : m_fields )
  {
    globalIndex const adjustment = field.globalOffset - field.rankOffset;

    arrayView1d< globalIndex > const indexArray =
      field.manager->getReference< array1d< globalIndex > >( field.key ).toView();

    forAll< parallelHostPolicy >( field.manager->numOwnedObjects(), [indexArray, adjustment]( localIndex const i )
    {
      indexArray[i] += adjustment;
    } );

    string_array fieldNames;
    fieldNames.emplace_back( field.key );
    CommunicationTools::getInstance().synchronizeFields( fieldNames, *field.manager, m_domain->getNeighbors(), false );
  }
}

template< typename MATRIX >
void DofManager::makeRestrictor( string const & fieldName,
                                 MPI_Comm const & comm,
                                 bool const transpose,
                                 MATRIX & restrictor ) const
{
  GEOS_ERROR_IF( !m_reordered, "Cannot make restrictors before reorderByRank() has been called." );

  FieldDescription const & field = getField( fieldName );
  restrictor.createWithLocalSize( transpose ? numLocalDofs() : field.numLocalDof,
                                  transpose ? field.numLocalDof : numLocalDofs(),
                                  1,
                                  comm );
  restrictor.open();

  array1d< globalIndex > rows( field.numLocalDof );
  array1d< globalIndex > cols( field.numLocalDof );

  forAll< parallelHostPolicy >( field.numLocalDof, [&field, transpose,
                                                    rows = rows.toView(),
                                                    cols = cols.toView()]( localIndex const i )
  {
    globalIndex row = field.rankOffset + i;
    globalIndex col = field.globalOffset + i;
    if( transpose )
    {
      std::swap( row, col );
    }
    rows[i] = row;
    cols[i] = col;
  } );

  array1d< real64 > vals( field.numLocalDof );
  vals.setValues< parallelHostPolicy >( 1.0 );

  restrictor.insert( rows.toViewConst(),
                     cols.toViewConst(),
                     vals.toViewConst() );

  restrictor.close();
}

#define MAKE_DOFMANAGER_METHOD_INST( LAI ) \
  template void DofManager::makeRestrictor( string const & fieldName, \
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

} // geosx
} // multiscale
