/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

namespace geosx
{

// .... DOF MANAGER :: CONSTRUCTOR
DofManager::DofManager()
{
  mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX );
  mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );

  // we pre-allocate an oversized array to store connectivity type
  // instead of resizing it dynamically as fields are added.
  m_connectivity.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );

  for( localIndex i = 0 ; i < MAX_NUM_FIELDS ; ++i )
  {
    for( localIndex j = 0 ; j < MAX_NUM_FIELDS ; ++j )
    {
      m_connectivity[i][j] = Connectivity::None;
    }
  }

  // we pre-allocate an oversized array to store sparsity pattern type
  // instead of resizing it dynamically as fields are added.
  m_sparsityPattern.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );
  for( localIndex i = 0 ; i < MAX_NUM_FIELDS ; ++i )
  {
    for( localIndex j = 0 ; j < MAX_NUM_FIELDS ; ++j )
    {
      m_sparsityPattern[i][j] = std::make_pair( nullptr, nullptr );
    }
  }
}

// .... DOF MANAGER :: SET MESH
void DofManager::setMesh( DomainPartition * const domain,
                          localIndex const meshLevelIndex,
                          localIndex const meshBodyIndex )
{
  GEOS_ERROR_IF( m_meshLevel != nullptr, "A mesh is already assigned to this DofManager." );
  m_domain = domain;
  m_meshLevel = m_domain->getMeshBodies()->GetGroup<MeshBody>( meshBodyIndex )->
                getMeshLevel( integer_conversion<int>( meshLevelIndex ) );
}

// .... DOF MANAGER :: FIELD INDEX
localIndex DofManager::fieldIndex( string const & key ) const
{
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
  {
    if( m_fields[i].name == key )
    {
      return i;
    }
  }
  GEOS_ERROR( "Field's string key not found in list of active fields." );
  return -1;
}

// .... DOF MANAGER :: KEY IN USE
bool DofManager::keyInUse( string const & key ) const
{
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
  {
    if( m_fields[i].name == key )
    {
      return true;
    }
  }
  return false;
}

// Return global number of dofs across all processors. If field argument is empty, return monolithic size.
globalIndex DofManager::numGlobalDofs( string const & field ) const
{
  if( field.length() > 0 )
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( field ), "numGlobalDofs: requested field name must be already existing." );

    // get field index
    localIndex fieldIdx = fieldIndex( field );

    return m_fields[fieldIdx].numGlobalRows;
  }
  else
  {
    globalIndex sumGlobalDofs = 0;
    for( localIndex i = 0 ; i < m_fields.size() ; ++i )
    {
      sumGlobalDofs += m_fields[i].numGlobalRows;
    }
    return sumGlobalDofs;
  }
}

// Return local number of dofs across all processors. If field argument is empty, return monolithic size.
localIndex DofManager::numLocalDofs( string const & field ) const
{
  if( field.length() > 0 )
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( field ), "numLocalDofs: requested field name must be already existing." );

    // get field index
    localIndex fieldIdx = fieldIndex( field );

    return m_fields[fieldIdx].numLocalRows;
  }
  else
  {
    localIndex sumLocalDofs = 0;
    for( localIndex i = 0 ; i < m_fields.size() ; ++i )
    {
      sumLocalDofs += m_fields[i].numLocalRows;
    }
    return sumLocalDofs;
  }
}

// Just an interface to allow only three parameters
void DofManager::addField( string const & field,
                           Location const location,
                           Connectivity const connectivity )
{
  addField( field, location, connectivity, 1, array1d<string>() );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addField( string const & field,
                           Location const location,
                           Connectivity const connectivity,
                           localIndex const components )
{
  addField( field, location, connectivity, components, array1d<string>() );
}

// Just another interface to allow four parameters (no components)
void DofManager::addField( string const & field,
                           Location const location,
                           Connectivity const connectivity,
                           string_array const & regions )
{
  addField( field, location, connectivity, 1, regions );
}

// The real function, allowing the creation of self-connected blocks
void DofManager::addField( string const & field,
                           Location const location,
                           Connectivity const connectivity,
                           localIndex const components,
                           string_array const & regions )
{
  // check if the field name is already being used
  GEOS_ERROR_IF( keyInUse( field ), "Requested field name matches an existing field in the DofManager." );

  // save field description to list of active fields
  FieldDescription description;

  description.name = field;
  description.location = location;
  description.numComponents = components;
  description.key = field + "_dof_indices";
  description.docstring = field + " dof indices";

  if( components > 1 )
  {
    description.docstring += " (with " + std::to_string( components ) + "-component blocks)";
  }

  // save pointers to "active" element regions
  ElementRegionManager * const elemManager = m_meshLevel->getElemManager();

  // retrieve full list of regions
  if( regions.size() == 0 )
  {
    auto const & regionListPtr = elemManager->GetRegions().keys();
    string_array regionNames( regionListPtr.size() );

    for( auto& regionPtr : regionListPtr )
    {
      regionNames[regionPtr.second] = regionPtr.first;
    }

    description.regionNames = regionNames;
  }
  else
  {
    description.regionNames = regions;
  }

  m_fields.push_back( description );

  localIndex numFields = m_fields.size();
  GEOS_ERROR_IF( numFields > MAX_NUM_FIELDS, "Limit on DofManager's MAX_NUM_FIELDS exceeded." );

  // temp reference to last field
  FieldDescription & last = m_fields[numFields - 1];

  localIndex numTotalRegions = elemManager->numRegions();
  localIndex numActiveRegions = regions.size() == 0 ? numTotalRegions : regions.size();

  last.regionPtrs.resize( numActiveRegions );
  for( localIndex er = 0 ; er < numActiveRegions ; ++er )
  {
    // Get region by name
    last.regionPtrs[er] = elemManager->GetRegion( last.regionNames[er] );
    GEOS_ERROR_IF( last.regionPtrs[er] == nullptr, "Specified element region not found" );
  }

  // based on location, allocate an index array for this field
  switch( location )
  {
    case Location::Elem:
      createIndexArray_ElemVersion( last );
      break;
    case Location::Face:
      createIndexArray_NodeOrFaceVersion( last );
      break;
    case Location::Node:
      createIndexArray_NodeOrFaceVersion( last );
      break;
    default:
      GEOS_ERROR( "DoF support location is not yet supported" );
  }

  // determine field's global offset
  if( numFields > 1 )
  {
    FieldDescription & prev = m_fields[numFields - 2];
    last.fieldOffset = prev.fieldOffset + prev.numGlobalRows;
  }
  else
  {
    last.fieldOffset = 0;
  }

  // save field's connectivity type (self-to-self)
  m_connectivity[numFields - 1][numFields - 1] = connectivity;

  // add sparsity pattern (LC matrix)
  Dof_SparsityPattern connLocPattLocal;
  addDiagSparsityPattern( connLocPattLocal, numFields - 1, connectivity );

  // TRILINOS interface
  localIndex maxEntriesPerRow = 0;
  for( localIndex i = 0 ; i < connLocPattLocal.nRows ; ++i )
  {
    maxEntriesPerRow = std::max( maxEntriesPerRow,
                                 connLocPattLocal.rowLengths[i + 1] - connLocPattLocal.rowLengths[i] );
  }

  last.connLocPattern = new ParallelMatrix();
  ParallelMatrix* connLocPattDistr = last.connLocPattern;
  connLocPattDistr->createWithGlobalSize( connLocPattLocal.nRows,
                                          connLocPattLocal.nCols,
                                          maxEntriesPerRow,
                                          MPI_COMM_GEOSX );
  for( globalIndex i = 0 ; i < connLocPattLocal.nRows ; ++i )
  {
    localIndex nnz = connLocPattLocal.rowLengths[i + 1] - connLocPattLocal.rowLengths[i];
    if( nnz > 0 )
    {
      real64_array values( nnz );
      values = 1;
      connLocPattDistr->insert( i,
                                connLocPattLocal.colIndices.data( connLocPattLocal.rowLengths[i] ),
                                values.data(),
                                nnz );
    }
  }
  connLocPattDistr->close();

  // log some basic info
  GEOS_LOG_RANK_0( "DofManager :: Added field .... " << last.docstring );
  GEOS_LOG_RANK_0( "DofManager :: Global dofs .... " << last.numGlobalRows );
  GEOS_LOG_RANK_0( "DofManager :: Field offset ... " << last.fieldOffset );
}

// addField: allow the usage of a predefine location-connection pattern (user-defined)
// Interface to allow only two parameters
void DofManager::addField( string const & field,
                           ParallelMatrix const & connLocInput )
{
  DofManager::addField( field,
                        connLocInput,
                        1,
                        Connectivity::USER_DEFINED );
}

// Just another interface to allow three parameters (no connectivity)
void DofManager::addField( string const & field,
                           ParallelMatrix const & connLocInput,
                           localIndex const components )
{
  DofManager::addField( field,
                        connLocInput,
                        components,
                        Connectivity::USER_DEFINED );
}

// Just another interface to allow three parameters (no components)
void DofManager::addField( string const & field,
                           ParallelMatrix const & connLocInput,
                           Connectivity const connectivity )
{
  DofManager::addField( field,
                        connLocInput,
                        1,
                        connectivity );
}

// The real function
void DofManager::addField( string const & field,
                           ParallelMatrix const & connLocInput,
                           localIndex const components,
                           Connectivity const connectivity )
{
  // check if the field name is already being used
  GEOS_ERROR_IF( keyInUse( field ), "Requested field name matches an existing field in the DofManager." );

  // save field description to list of active fields
  FieldDescription description;

  description.name = field;
  description.location = Location::USER_DEFINED;
  description.numComponents = components;
  description.key = field + "_dof_indices";
  description.docstring = field + " dof indices";

  if( components > 1 )
  {
    description.docstring += " (with " + std::to_string( components ) + "-component blocks)";
  }

  m_fields.push_back( description );

  localIndex numFields = m_fields.size();
  GEOS_ERROR_IF( numFields > MAX_NUM_FIELDS, "Limit on DofManager's MAX_NUM_FIELDS exceeded." );

  // temp reference to last field
  FieldDescription & last = m_fields[numFields - 1];

  // determine field's global offset
  if( numFields > 1 )
  {
    FieldDescription & prev = m_fields[numFields - 2];
    last.fieldOffset = prev.fieldOffset + prev.numGlobalRows;
  }
  else
  {
    last.fieldOffset = 0;
  }

  // save field's connectivity type (self-to-self)
  m_connectivity[numFields - 1][numFields - 1] = connectivity;

  // save the user-provided location-connectivity matrix
  last.connLocPattern = new ParallelMatrix(connLocInput);

  // compute useful values (number of local and global rows)
  last.numLocalRows = numMyCols(connLocInput);
  last.numLocalNodes = last.numLocalRows / components;

  localIndex_array localGather;

  CommunicationTools::allGather( last.numLocalNodes, localGather );

  last.numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    last.numGlobalRows += localGather[p];
  }
  last.numGlobalRows *= components;

  last.firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    last.firstLocalRow += localGather[p];
  }
  last.firstLocalRow *= components;

  if ( connectivity == Connectivity::Elem )
  {
    globalIndex_array globalGather;

    last.firstLocalConnectivity = connLocInput.unwrappedPointer()->NumMyRows();

    CommunicationTools::allGather( last.firstLocalConnectivity, globalGather );

    last.firstLocalConnectivity = 0;
    for( localIndex p = 0 ; p < mpiRank ; ++p )
    {
      last.firstLocalConnectivity += globalGather[p];
    }
  }
  else
  {
    last.firstLocalConnectivity = 0;
  }

  // log some basic info
  GEOS_LOG_RANK_0( "DofManager :: Added field .... " << last.docstring );
  GEOS_LOG_RANK_0( "DofManager :: Global dofs .... " << last.numGlobalRows );
  GEOS_LOG_RANK_0( "DofManager :: Field offset ... " << last.fieldOffset );
}

// .... DOF MANAGER :: CREATE INDEX ARRAY
void DofManager::createIndexArray_NodeOrFaceVersion( FieldDescription & field,
                                                     localIndex_array const & activeRegionsInput ) const
{
  // step 0. register an index array with default = LocationStatus::notAssigned
  ObjectManagerBase *
  baseManager = field.location == Location::Node ?
                static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() ) :
                static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() ) ;

  baseManager->RegisterViewWrapper<globalIndex_array>( field.key )->
    setApplyDefaultValue( static_cast<globalIndex>( LocationStatus::notAssigned ) )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )->
    setDescription( field.docstring );

  globalIndex_array & indexArray = baseManager->getReference<globalIndex_array>( field.key );

  // compute activeRegions
  localIndex_array activeRegions( activeRegionsInput );
  if( activeRegions.size() == 0 )
  {
    activeRegions.resize( field.regionPtrs.size() );
    activeRegions = 1;
  }

  // step 1. loop over all active regions
  //         determine number of local rows
  //         and sequentially number objects
  field.numLocalNodes = 0;
  localIndex numLocalNodesWithGhost = 0;
  globalIndex numLocalConnectivity = 0;

  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
  {
    if( activeRegions[er] >= 0 )
    {
      for( localIndex esr = 0 ; esr < field.regionPtrs[er]->numSubRegions() ; esr++ )
      {
        CellElementSubRegion const * const subRegion = field.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
        integer_array const & ghostRank = subRegion->m_ghostRank;

        localIndex_array2d const &
        map = field.location == Location::Node ?
              subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().nodeListString )->reference() :
              subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference() ;

        // Set which process owns the boundary nodes/faces
        for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
        {
          if( !( ghostRank[e] < 0 ) and ghostRank[e] < mpiRank )
          {
            for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
            {
              indexArray[map[e][n]] = static_cast<globalIndex>( LocationStatus::notMyGhostLocation );
            }
          }
        }

        for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
        {
          for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
          {
            localIndex i = map[e][n];
            if( indexArray[i] == static_cast<globalIndex>( LocationStatus::notAssigned ) )
            {
              indexArray[i] = numLocalNodesWithGhost;
              numLocalNodesWithGhost++;
              if( ghostRank[e] < 0 )
              {
                ++field.numLocalNodes;
              }
            }
          }
          ++numLocalConnectivity;
        }
      }
    }
  }

  // step 2. gather row counts across ranks
  localIndex_array localGather;

  CommunicationTools::allGather( field.numLocalNodes, localGather );

  field.numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    field.numGlobalRows += localGather[p];
  }

  field.firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    field.firstLocalRow += localGather[p];
  }

  // for starting the connectivity numbering
  globalIndex_array globalGather;
  CommunicationTools::allGather( numLocalConnectivity, globalGather );

  field.firstLocalConnectivity = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    field.firstLocalConnectivity += globalGather[p];
  }

  // step 3. adjust local values to reflect processor offset
  for( localIndex n = 0 ; n < indexArray.size() ; ++n )
  {
    if( indexArray[n] != static_cast<globalIndex>( LocationStatus::notAssigned ) )
    {
      indexArray[n] += field.firstLocalRow;
    }
  }

  // step 4. synchronize across ranks
  std::map<string, string_array> fieldNames;

  if( field.location == Location::Node )
  {
    fieldNames["node"].push_back( field.key );
  }
  else
  {
    fieldNames["face"].push_back( field.key );
  }

  CommunicationTools::
  SynchronizeFields( fieldNames, m_meshLevel,
                     m_domain->getReference<array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components
  field.numLocalRows = field.numLocalNodes * field.numComponents;
  field.numGlobalRows *= field.numComponents;
  field.firstLocalRow *= field.numComponents;

  // Replace LocationStatus::notMyGhostLocation with LocationStatus::notAssigned (if any)
  for( localIndex i = 0 ; i < indexArray.size() ; ++i )
  {
    if( indexArray[i] == static_cast<globalIndex>( LocationStatus::notMyGhostLocation ) )
    {
      indexArray[i] = static_cast<globalIndex>( LocationStatus::notAssigned );
    }
  }
}

// .... DOF MANAGER :: CREATE INDEX ARRAY :: ELEMENT VERSION
//      TODO: revise to look more like node version.
//            may even be able to condense to one function.
void DofManager::createIndexArray_ElemVersion( FieldDescription & field ) const
{
  // step 1. loop over all active regions
  //         determine number of local rows
  field.numLocalNodes = 0;
  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
  {
    field.regionPtrs[er]->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    {
      localIndex numGhost = subRegion->GetNumberOfGhosts();
      field.numLocalNodes += subRegion->size() - numGhost;
    } );
  }

  // step 2. gather row counts across ranks
  localIndex_array localGather;

  CommunicationTools::allGather( field.numLocalNodes, localGather );

  field.numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    field.numGlobalRows += localGather[p];
  }

  field.firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    field.firstLocalRow += localGather[p];
  }

  // step 3. loop again (sequential policy)
  //         allocate the index array
  //         set unique global indices
  globalIndex count = 0;

  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
  {
    for( localIndex esr = 0 ; esr < field.regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellElementSubRegion * const subRegion = field.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );

      subRegion->RegisterViewWrapper<globalIndex_array>( field.key )->
        setApplyDefaultValue( static_cast<globalIndex>( LocationStatus::notAssigned ) )->
        setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )->
        setDescription( field.docstring );

      globalIndex_array & indexArray = subRegion->getReference<globalIndex_array>( field.key );
      integer_array const & ghostRank = subRegion->m_ghostRank;

      GEOS_ERROR_IF( indexArray.size() != ghostRank.size(), "Mismatch in ghost rank and index array sizes." );

      for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
      {
        if( ghostRank[elem] < 0 )
        {
          indexArray[elem] = field.firstLocalRow + count;
          count++;
        }
      }
    }
  }

  GEOS_ERROR_IF( count != field.numLocalNodes, "Mismatch during assignment of local row indices" );

  // for starting the connectivity numbering
  globalIndex_array globalGather;

  CommunicationTools::allGather( count, globalGather );

  field.firstLocalConnectivity = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    field.firstLocalConnectivity += globalGather[p];
  }

  // step 4. synchronize across ranks
  std::map<string, string_array> fieldNames;
  fieldNames["elem"].push_back( field.key );

  CommunicationTools::
  SynchronizeFields( fieldNames, m_meshLevel,
                     m_domain->getReference<array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components
  field.numLocalRows = field.numLocalNodes * field.numComponents;
  field.numGlobalRows *= field.numComponents;
  field.firstLocalRow *= field.numComponents;
}

// Create the sparsity pattern (location-location). High level interface
void DofManager::getSparsityPattern( ParallelMatrix & locLocDistr,
                                     string const & rowField,
                                     string const & colField ) const
{
  localIndex rowFieldIndex, colFieldIndex;

  if( rowField.length() > 0 )
  {
    // check if the row field name is already added
    GEOS_ERROR_IF( !keyInUse( rowField ), "getSparsityPattern: requested field name must be already existing." );

    // get row field index
    rowFieldIndex = fieldIndex( rowField );
  }
  else
  {
    rowFieldIndex = -1;
  }

  if( colField.length() > 0 )
  {
    // check if the col field name is already added
    GEOS_ERROR_IF( !keyInUse( colField ), "getSparsityPattern: requested field name must be already existing." );

    // get col field index
    colFieldIndex = fieldIndex( colField );
  }
  else
  {
    colFieldIndex = -1;
  }

  if( rowFieldIndex * colFieldIndex < 0 )
  {
    GEOS_ERROR( "getSparsityPattern accepts both two field names and none, instead just one is provided." );
  }

  // Call the low level routine
  getSparsityPattern( locLocDistr, rowFieldIndex, colFieldIndex );
}

// Create the sparsity pattern (location-location). Low level interface
void DofManager::getSparsityPattern( ParallelMatrix & locLocDistr,
                                     localIndex const rowFieldIndex,
                                     localIndex const colFieldIndex ) const
{
  GEOS_ERROR_IF( rowFieldIndex * colFieldIndex < 0,
                 "getSparsityPattern accepts both two existing field indices (positive values) and "
                 "two negative values (entire Jacobian rows/columns), instead just one index is positive.");

  if( rowFieldIndex >= 0 and colFieldIndex == rowFieldIndex )
  {
    // Diagonal block
    ParallelMatrix const * connLocPattDistr = m_fields[rowFieldIndex].connLocPattern;

    locLocDistr.createWithGlobalSize( connLocPattDistr->globalCols(), 1, MPI_COMM_GEOSX );
    MatrixMatrixMultiply( *connLocPattDistr, true, *connLocPattDistr, false, locLocDistr );
  }
  else if( rowFieldIndex >= 0 and colFieldIndex >= 0 )
  {
    // ExtraDiagonal (coupling) block
    if( m_connectivity[rowFieldIndex][colFieldIndex] != Connectivity::None )
    {
      if( m_sparsityPattern[rowFieldIndex][colFieldIndex].first != nullptr )
      {
        ParallelMatrix const * rowConnLocPattDistr = m_sparsityPattern[rowFieldIndex][colFieldIndex].first;
        ParallelMatrix const * colConnLocPattDistr = m_sparsityPattern[rowFieldIndex][colFieldIndex].second;

        locLocDistr.createWithGlobalSize( rowConnLocPattDistr->globalCols(),
                                          colConnLocPattDistr->globalCols(),
                                          1,
                                          MPI_COMM_GEOSX );
        MatrixMatrixMultiply( *rowConnLocPattDistr, true, *colConnLocPattDistr, false, locLocDistr, false );
      }
      else
      {
        ParallelMatrix const * rowConnLocPattDistr = m_sparsityPattern[colFieldIndex][rowFieldIndex].first;
        ParallelMatrix const * colConnLocPattDistr = m_sparsityPattern[colFieldIndex][rowFieldIndex].second;

        locLocDistr.createWithGlobalSize( colConnLocPattDistr->globalCols(),
                                          rowConnLocPattDistr->globalCols(),
                                          1,
                                          MPI_COMM_GEOSX );
        MatrixMatrixMultiply( *colConnLocPattDistr, true, *rowConnLocPattDistr, false, locLocDistr, false );
      }
    }
    else
    {
      // Empty block
      globalIndex nRows = m_fields[rowFieldIndex].connLocPattern->globalCols();
      globalIndex nCols = m_fields[colFieldIndex].connLocPattern->globalCols();
      locLocDistr.createWithGlobalSize( nRows, nCols, 0, MPI_COMM_GEOSX );
    }
  }
  else if( rowFieldIndex < 0 and colFieldIndex < 0 )
  {
    // Create the global matrix
    globalIndex sumGlobalDofs = 0;
    for( localIndex i = 0 ; i < m_fields.size() ; ++i )
    {
      sumGlobalDofs += m_fields[i].numGlobalRows;
    }
    locLocDistr.createWithGlobalSize( sumGlobalDofs, sumGlobalDofs, 1, MPI_COMM_GEOSX );

    ParallelMatrix localPattern;

    // Loop over all fields
    for( localIndex iGlo = 0 ; iGlo < m_fields.size() ; ++iGlo )
    {
      // Loop over all fields
      for( localIndex jGlo = 0 ; jGlo < m_fields.size() ; ++jGlo )
      {
        if( iGlo == jGlo )
        {
          // Diagonal block
          ParallelMatrix const * connLocPattDistr = m_fields[iGlo].connLocPattern;

          localPattern.createWithGlobalSize( connLocPattDistr->globalCols(), 1, MPI_COMM_GEOSX );
          MatrixMatrixMultiply( *connLocPattDistr, true, *connLocPattDistr, false, localPattern );

          for( globalIndex i = localPattern.ilower() ; i < localPattern.iupper() ; ++i )
          {
            globalIndex_array indices;
            real64_array values;
            localPattern.getRowCopy( i, indices, values );
            if( indices.size() > 0 )
            {
              for( globalIndex j = 0 ; j < indices.size() ; ++j )
              {
                indices[j] += m_fields[jGlo].fieldOffset;
              }
              locLocDistr.insert( i + m_fields[iGlo].fieldOffset, indices, values );
            }
          }
        }
        else if( m_connectivity[iGlo][jGlo] != Connectivity::None )
        {
          // ExtraDiagonal (coupling) block
          if( m_sparsityPattern[iGlo][jGlo].first != nullptr )
          {
            ParallelMatrix const * rowConnLocPattDistr = m_sparsityPattern[iGlo][jGlo].first;
            ParallelMatrix const * colConnLocPattDistr = m_sparsityPattern[iGlo][jGlo].second;

            localPattern.createWithGlobalSize( rowConnLocPattDistr->globalCols(), colConnLocPattDistr->globalCols(),
                                               1, MPI_COMM_GEOSX );
            MatrixMatrixMultiply( *rowConnLocPattDistr, true, *colConnLocPattDistr, false, localPattern, false );
          }
          else
          {
            ParallelMatrix const * rowConnLocPattDistr = m_sparsityPattern[jGlo][iGlo].first;
            ParallelMatrix const * colConnLocPattDistr = m_sparsityPattern[jGlo][iGlo].second;

            localPattern.createWithGlobalSize( colConnLocPattDistr->globalCols(), rowConnLocPattDistr->globalCols(),
                                               1, MPI_COMM_GEOSX );
            MatrixMatrixMultiply( *colConnLocPattDistr, true, *rowConnLocPattDistr, false, localPattern, false );
          }

          // Assembly (with right offsets)
          for( globalIndex i = localPattern.ilower() ; i < localPattern.iupper() ; ++i )
          {
            globalIndex_array indices;
            real64_array values;
            localPattern.getRowCopy( i, indices, values );
            if( indices.size() > 0 )
            {
              for( globalIndex j = 0 ; j < indices.size() ; ++j )
              {
                indices[j] += m_fields[jGlo].fieldOffset;
              }
              locLocDistr.insert( i + m_fields[iGlo].fieldOffset, indices, values );
            }
          }
        }
      }
    }
  }
  locLocDistr.close();
}

// Permute the GLOBAL sparsity pattern (location-location). Low level interface
void DofManager::permuteSparsityPattern( ParallelMatrix const & locLocDistr,
                                         ParallelMatrix const & permutation,
                                         ParallelMatrix & permutedMatrix ) const
{
  // Performe the product B = P^t*A*P
  // Matrix C = A*P
  ParallelMatrix productStep1;

  productStep1.createWithGlobalSize( locLocDistr.globalRows(), permutation.globalCols(), 1, MPI_COMM_GEOSX );
  MatrixMatrixMultiply( locLocDistr, false, permutation, false, productStep1 );

  // Matrix B = P*C
  permutedMatrix.createWithGlobalSize( permutation.globalCols(), productStep1.globalCols(), 1, MPI_COMM_GEOSX );
  MatrixMatrixMultiply( permutation, true, productStep1, false, permutedMatrix );
}

// Just an interface to allow only three parameters
void DofManager::addCoupling( string const & rowField,
                              string const & colField,
                              Connectivity const connectivity ) const
{
  addCoupling( rowField, colField, connectivity, string_array(), true );
}

// Just another interface to allow four parameters (no symmetry)
void DofManager::addCoupling( string const & rowField,
                              string const & colField,
                              Connectivity const connectivity,
                              string_array const & regions = string_array() ) const
{
  addCoupling( rowField, colField, connectivity, regions, true );
}

// Just another interface to allow four parameters (no regions)
void DofManager::addCoupling( string const & rowField,
                              string const & colField,
                              Connectivity const connectivity,
                              bool const symmetric ) const
{
  addCoupling( rowField, colField, connectivity, string_array(), symmetric );
}

// The real function, allowing the creation of coupling blocks
void DofManager::addCoupling( string const & rowField,
                              string const & colField,
                              Connectivity const connectivity,
                              string_array const & regions,
                              bool const symmetric ) const
{
  // check if the row field name is already added
  GEOS_ERROR_IF( !keyInUse( rowField ), "addCoupling: requested field name must be already existing." );

  // check if the col field name is already added
  GEOS_ERROR_IF( !keyInUse( colField ), "addCoupling: requested field name must be already existing." );

  // get row field index
  localIndex const rowFieldIndex = fieldIndex( rowField );

  // get col field index
  localIndex const colFieldIndex = fieldIndex( colField );

  // Row field description
  FieldDescription const & rowFieldDesc = m_fields[rowFieldIndex];

  // Col field description
  FieldDescription const & colFieldDesc = m_fields[colFieldIndex];

  // Row field regionName
  array1d<string> const rowFieldRegionNames = rowFieldDesc.regionNames;

  // Col field regionName
  array1d<string> const colFieldRegionNames = colFieldDesc.regionNames;

  string_array regionsList( regions );
  if( regionsList.size() == 0 )
  {
    // Detect common regions between row and col fields
    // Sort list of regions (for row and col fields)
    string_array rowFieldRegionNamesSorted( rowFieldRegionNames );
    string_array colFieldRegionNamesSorted( colFieldRegionNames );

    std::sort( rowFieldRegionNamesSorted.begin(), rowFieldRegionNamesSorted.end() );
    std::sort( colFieldRegionNamesSorted.begin(), colFieldRegionNamesSorted.end() );

    // Resize regionsList
    regionsList.resize( rowFieldRegionNames.size() + colFieldRegionNames.size() );

    // Find common regions
    string_array::iterator
    it = std::set_intersection( rowFieldRegionNamesSorted.begin(), rowFieldRegionNamesSorted.end(),
                                colFieldRegionNamesSorted.begin(), colFieldRegionNamesSorted.end(), regionsList.begin() );
    regionsList.resize( it - regionsList.begin() );
  }

  bool areDefinedRegions = false;
  localIndex_array rowFieldRegionIndex;
  localIndex_array colFieldRegionIndex;
  if( regionsList.size() <= std::min( rowFieldRegionNames.size(), colFieldRegionNames.size() ) )
  {
    areDefinedRegions = true;
    rowFieldRegionIndex.resize( rowFieldRegionNames.size() );
    colFieldRegionIndex.resize( colFieldRegionNames.size() );
    rowFieldRegionIndex = -1;
    colFieldRegionIndex = -1;
    // Check if row and col regions are the same
    localIndex rowNameIndex = 0;
    localIndex colNameIndex = 0;

    for( array1d<string>::const_iterator regionName = regionsList.begin() ; regionName != regionsList.end() ;
        ++regionName )
    {
      localIndex
      rowID = std::find( rowFieldRegionNames.begin(), rowFieldRegionNames.end(), *regionName )
            - rowFieldRegionNames.begin();
      bool rowDefined = rowID < rowFieldRegionNames.size();
      if( rowDefined )
      {
        rowFieldRegionIndex[rowID] = rowNameIndex++;
      }
      areDefinedRegions &= rowDefined;
      localIndex
      colID = std::find( colFieldRegionNames.begin(), colFieldRegionNames.end(), *regionName )
            - colFieldRegionNames.begin();
      bool colDefined = colID < colFieldRegionNames.size();
      if( colDefined )
      {
        colFieldRegionIndex[colID] = colNameIndex++;
      }
      areDefinedRegions &= colDefined;
    }

    if( !areDefinedRegions )
    {
      GEOS_ERROR_IF( !areDefinedRegions,
                     "addCoupling: regions where coupling is defined have to belong to already existing fields." );
      return;
    }
  }
  else
  {
    GEOS_ERROR_IF( !areDefinedRegions,
                   "addCoupling: regions where coupling is defined have to belong to already existing fields." );
    return;
  }

  // Check if already defined
  if( m_connectivity[rowFieldIndex][colFieldIndex] != Connectivity::None )
  {
    GEOS_ERROR( "addCoupling: coupling already defined with another connectivity." );
    return;
  }

  // save field's connectivity type (rowField to colField)
  m_connectivity[rowFieldIndex][colFieldIndex] = connectivity;

  if( connectivity != Connectivity::None )
  {
    // get pointer to the right matrix pair
    ParallelMatrix *& rowPattern = m_sparsityPattern[rowFieldIndex][colFieldIndex].first;
    ParallelMatrix *& colPattern = m_sparsityPattern[rowFieldIndex][colFieldIndex].second;

    // add sparsity pattern
    addExtraDiagSparsityPattern( rowPattern, colPattern, rowFieldIndex, colFieldIndex, rowFieldRegionIndex,
                                 colFieldRegionIndex, connectivity );
  }

  // Set connectivity with active symmetry flag
  if( symmetric )
  {
    m_connectivity[colFieldIndex][rowFieldIndex] = connectivity;
  }
}

// Get global indices for dofs connected by the connectivity type in the given
// region/subregion combination
void DofManager::getIndices( globalIndex_array & indices,
                             Connectivity const connectivity,
                             localIndex const region,
                             localIndex const subregion,
                             localIndex const index,
                             string const & field ) const
{
  // resize to 0
  indices.resize( 0 );

  if( connectivity != Connectivity::Elem )
  {
    getIndices( indices, connectivity, index, field );
  }
  else
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( field ), "getIndices: requested field name must be already existing." );

    // get field index
    localIndex const fieldIdx = fieldIndex( field );

    // check connectivity
    if( m_connectivity[fieldIdx][fieldIdx] == connectivity )
    {
      // get field description
      FieldDescription const & fieldDesc = m_fields[fieldIdx];

      globalIndex_array firstLocalConnectivity( fieldDesc.regionPtrs.size() );

      if( region < fieldDesc.regionPtrs.size() )
      {
        if( subregion < fieldDesc.regionPtrs[region]->numSubRegions() )
        {
          for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
          {
            localIndex localCount = 0;
            for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
            {
              CellElementSubRegion const * const
              subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
              integer_array const & ghostRank = subRegion->m_ghostRank;

              for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
              {
                if( ghostRank[elem] < 0 )
                {
                  ++localCount;
                }
              }
            }

            localIndex_array localGather;
            CommunicationTools::allGather( localCount, localGather );

            firstLocalConnectivity[er] = 0;
            for( localIndex p = 0 ; p < mpiRank ; ++p )
            {
              firstLocalConnectivity[er] += localGather[p];
            }
          }
        }
      }

      if( region < fieldDesc.regionPtrs.size() )
      {
        if( subregion < fieldDesc.regionPtrs[region]->numSubRegions() )
        {
          globalIndex globalCount = fieldDesc.firstLocalConnectivity;
          for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
          {
            globalIndex localCount = firstLocalConnectivity[er];
            for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
            {
              CellElementSubRegion const * const
              subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
              integer_array const & ghostRank = subRegion->m_ghostRank;

              for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
              {
                if( ghostRank[elem] < 0 )
                {
                  if( er == region && esr == subregion && index == localCount )
                  {
                    localIndex localRow = ParallelMatrixGetLocalRowID( *fieldDesc.connLocPattern, globalCount );
                    if( localRow >= 0 )
                    {
                      // Retrieve row
                      real64_array values;
                      fieldDesc.connLocPattern->getRowCopy( globalCount, indices, values );
                      // Add offset
                      for( localIndex i = 0 ; i < indices.size() ; ++i )
                      {
                        indices[i] += m_fields[fieldIdx].fieldOffset;
                      }
                    }
                  }
                  ++localCount;
                  ++globalCount;
                }
              }
            }
          }
        }
      }
    }
  }
}

// Get global indices for dofs connected by the connectivity type.
void DofManager::getIndices( globalIndex_array & indices,
                             Connectivity const connectivity,
                             localIndex const index,
                             string const & field ) const
{
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( field ), "getIndices: requested field name must be already existing." );

  // get field index
  const localIndex fieldIdx = fieldIndex( field );

  // resize to 0
  indices.resize( 0 );

  // check connectivity
  if( m_connectivity[fieldIdx][fieldIdx] == connectivity )
  {
    // get field description
    const FieldDescription & fieldDesc = m_fields[fieldIdx];

    if( index >= fieldDesc.connLocPattern->ilower() and index < fieldDesc.connLocPattern->iupper() )
    {
      // Retrieve row
      real64_array values;
      fieldDesc.connLocPattern->getRowCopy( index, indices, values );
      // Add offset
      for( localIndex i = 0 ; i < indices.size() ; ++i )
      {
        indices[i] += m_fields[fieldIdx].fieldOffset;
      }
    }
  }
}

// Create the permutation that collects together all DoFs of each MPI process
void DofManager::createPermutation( ParallelMatrix & permutation ) const
{
  globalIndex sumGlobalDofs = 0;
  globalIndex_array ilower( m_fields.size() );
  globalIndex_array iupper( m_fields.size() );
  globalIndex offset = 0;
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
  {
    sumGlobalDofs += m_fields[i].numGlobalRows;
    ParallelMatrix pattern;
    getSparsityPattern( pattern, i, i );
    ilower[i] = pattern.ilower();
    iupper[i] = pattern.iupper();
    offset += iupper[i] - ilower[i];
  }

  // Create parallel vector with global size
  permutation.createWithGlobalSize( sumGlobalDofs, sumGlobalDofs, 1, MPI_COMM_GEOSX );

  globalIndex_array globalGather;
  CommunicationTools::allGather( offset, globalGather );

  offset = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    offset += globalGather[p];
  }

  globalIndex k = offset;
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
  {
    globalIndex fieldOffset = m_fields[i].fieldOffset;
    for( globalIndex j = ilower[i] ; j < iupper[i] ; ++j )
    {
      permutation.insert( fieldOffset + j, k++, 1 );
    }
  }
  permutation.close();
}

// Add the specified sparsity pattern between rowField and colField to m_sparsityPattern
// collection
void DofManager::addExtraDiagSparsityPattern( ParallelMatrix *& rowConnLocPattDistr,
                                              ParallelMatrix *& colConnLocPattDistr,
                                              localIndex const & rowFieldIndex,
                                              localIndex const & colFieldIndex,
                                              localIndex_array const & rowActiveRegions,
                                              localIndex_array const & colActiveRegions,
                                              Connectivity const connectivity ) const
{
  // Row field description
  FieldDescription const & rowFieldDesc = m_fields[rowFieldIndex];

  // Col field description
  FieldDescription const & colFieldDesc = m_fields[colFieldIndex];

  // Compute the CL matrices for row and col fields
  Dof_SparsityPattern rowPatternLocal, colPatternLocal;
  addDiagSparsityPattern( rowPatternLocal, rowFieldIndex, connectivity, rowActiveRegions );
  addDiagSparsityPattern( colPatternLocal, colFieldIndex, connectivity, colActiveRegions );

  // TRILINOS interface

  // Pattern of connections and row locations
  localIndex maxEntriesPerRow = 0;
  for( localIndex i = 0 ; i < rowPatternLocal.nRows ; ++i )
  {
    maxEntriesPerRow = std::max( maxEntriesPerRow,
                                 rowPatternLocal.rowLengths[i + 1] - rowPatternLocal.rowLengths[i] );
  }

  rowConnLocPattDistr = new ParallelMatrix();
  rowConnLocPattDistr->createWithGlobalSize( rowPatternLocal.nRows,
                                             rowPatternLocal.nCols,
                                             maxEntriesPerRow,
                                             MPI_COMM_GEOSX );
  for( globalIndex i = 0 ; i < rowPatternLocal.nRows ; ++i )
  {
    localIndex nnz = rowPatternLocal.rowLengths[i + 1] - rowPatternLocal.rowLengths[i];
    if( nnz > 0 )
    {
      real64_array values( nnz );
      values = 1;
      rowConnLocPattDistr->insert( i,
                                   rowPatternLocal.colIndices.data( rowPatternLocal.rowLengths[i] ),
                                   values.data(),
                                   nnz );
    }
  }
  rowConnLocPattDistr->close();

  // Pattern of connections and column locations
  colConnLocPattDistr = new ParallelMatrix();
  colConnLocPattDistr->createWithGlobalSize( colPatternLocal.nRows,
                                             colPatternLocal.nCols,
                                             maxEntriesPerRow,
                                             MPI_COMM_GEOSX );
  for( globalIndex i = 0 ; i < colPatternLocal.nRows ; ++i )
  {
    localIndex nnz = colPatternLocal.rowLengths[i + 1] - colPatternLocal.rowLengths[i];
    if( nnz > 0 )
    {
      real64_array values( nnz );
      values = 1;
      colConnLocPattDistr->insert( i,
                                   colPatternLocal.colIndices.data( colPatternLocal.rowLengths[i] ),
                                   values.data(),
                                   nnz );
    }
  }
  colConnLocPattDistr->close();
}

// Compute the sparsity pattern of the matrix connectivity - location for the specified
// field (diagonal entry in the m_connectivity collection)
void DofManager::addDiagSparsityPattern( Dof_SparsityPattern & connLocPatt,
                                         localIndex const & fieldIdx,
                                         Connectivity const connectivity,
                                         localIndex_array const & activeRegionsInput ) const
{
  // get field description
  FieldDescription const & fieldDesc = m_fields[fieldIdx];

  // array to store the matrix in COO format
  array1d<indexPair> pairs;

  // detect active regions
  localIndex_array activeRegions( activeRegionsInput );
  if( activeRegionsInput.size() == 0 )
  {
    activeRegions.resize( fieldDesc.regionPtrs.size() );
    activeRegions = 1;
  }

  if( connectivity == Connectivity::None )
  {
    // Case of no connectivity (self connection)
    Connectivity selfConnectivity = static_cast<Connectivity>( fieldDesc.location );

    addDiagSparsityPattern( connLocPatt, fieldIdx, selfConnectivity, activeRegions );
    return;
  }
  else if( connectivity == Connectivity::Elem )
  {
    if( fieldDesc.location == Location::Elem )
    {
      // Case of connectivity = Elem and location = Elem
      globalIndex elemIndex = 0;
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            for( localIndex e = 0 ; e < ghostRank.size() ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                ++elemIndex;
              }
            }
          }
        }
      }

      globalIndex_array globalGather;
      CommunicationTools::allGather( elemIndex, globalGather );

      globalIndex firstLocalConnectivity = 0;
      for( localIndex p = 0 ; p < mpiRank ; ++p )
      {
        firstLocalConnectivity += globalGather[p];
      }

      elemIndex = 0;
      globalIndex count = 0;
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellElementSubRegion const * const
          subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;

          for( localIndex e = 0 ; e < ghostRank.size() ; ++e )
          {
            if( ghostRank[e] < 0 )
            {
              if( activeRegions[er] >= 0 )
              {
                for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                {
                  pairs.push_back( std::make_pair( firstLocalConnectivity + count,
                                                   fieldDesc.numComponents * ( fieldDesc.firstLocalConnectivity + elemIndex ) + i ) );
                }
                ++count;
              }
              ++elemIndex;
            }
          }
        }
      }
    }
    else
    {
      // Case of connectivity = Elem and location = Node or Face
      ObjectManagerBase *
      baseManager = fieldDesc.location == Location::Node ?
                    static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() ) :
                    static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() ) ;

      globalIndex_array & indexArray = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      globalIndex elemIndex = 0;
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            for( localIndex e = 0 ; e < ghostRank.size() ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                ++elemIndex;
              }
            }
          }
        }
      }

      globalIndex_array globalGather;
      CommunicationTools::allGather( elemIndex, globalGather );

      globalIndex firstLocalConnectivity = 0;
      for( localIndex p = 0 ; p < mpiRank ; ++p )
      {
        firstLocalConnectivity += globalGather[p];
      }

      elemIndex = 0;
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellElementSubRegion const * const
          subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;

          localIndex_array2d const &
          map = fieldDesc.location == Location::Node ?
                subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().nodeListString )->reference() :
                subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference() ;

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
          {
            if( ghostRank[e] < 0 )
            {
              if( activeRegions[er] >= 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  {
                    pairs.push_back( std::make_pair( firstLocalConnectivity + elemIndex,
                                                     fieldDesc.numComponents * indexArray[map[e][n]] + i ) );
                  }
                }
                ++elemIndex;
              }
            }
          }
        }
      }
    }
  }
  else if( connectivity == Connectivity::Face )
  {
    if( fieldDesc.location == Location::Face )
    {
      // Case of connectivity = Face and location = Face
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Face;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      FieldDescription fieldTmp2;
      fieldTmp2.regionNames = fieldDesc.regionNames;
      fieldTmp2.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp2.location = Location::Face;
      fieldTmp2.numComponents = 1;
      fieldTmp2.regionNames = fieldDesc.regionNames;
      fieldTmp2.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp2 );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFaceOrig = baseManager->getReference<globalIndex_array>( fieldTmp2.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Face
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference();

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  {
                    if( indexArrayFace[map[e][n]] >= 0 )
                    {
                      pairs.push_back( std::make_pair( indexArrayFace[map[e][n]],
                                                       fieldDesc.numComponents * indexArrayFaceOrig[map[e][n]] + i ) );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else if( fieldDesc.location == Location::Elem )
    {
      // Case of connectivity = Face and location = Elem
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Face;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Face
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;
            globalIndex_array const & indexArrayElem = subRegion->getReference<globalIndex_array>( fieldDesc.key );

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference();

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  {
                    if( indexArrayFace[map[e][n]] >= 0 )
                    {
                      pairs.push_back( std::make_pair( indexArrayFace[map[e][n]],
                                                       fieldDesc.numComponents * indexArrayElem[e] + i ) );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else
    {
      // Case of connectivity = Face and location = Node
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Face;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      globalIndex_array indexArrayFaceMarker( fieldTmp.numGlobalRows );
      indexArrayFaceMarker = 1;

      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference();
            localIndex_array nodeIndices;

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  globalIndex faceId = indexArrayFace[map[e][n]];
                  subRegion->GetFaceNodes( e, n, nodeIndices );
                  if( indexArrayFaceMarker( faceId ) >= 0 )
                  {
                    // Mark this face as already visited
                    indexArrayFaceMarker( faceId ) = -1;
                    for( localIndex j = 0 ; j < nodeIndices.size() ; ++j )
                    {
                      for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                      {
                        pairs.push_back( std::make_pair( faceId,
                                                         fieldDesc.numComponents * indexArrayNode[nodeIndices[j]] + i ) );
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
    if( fieldDesc.location == Location::Node )
    {
      // Case of connectivity = Node and location = Node
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Node;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      FieldDescription fieldTmp2;
      fieldTmp2.regionNames = fieldDesc.regionNames;
      fieldTmp2.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp2.location = Location::Node;
      fieldTmp2.numComponents = 1;
      fieldTmp2.regionNames = fieldDesc.regionNames;
      fieldTmp2.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp2 );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNodeOrig = baseManager->getReference<globalIndex_array>( fieldTmp2.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Node
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().nodeListString )->reference();

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  {
                    pairs.push_back( std::make_pair( indexArrayNode[map[e][n]],
                                                     fieldDesc.numComponents * indexArrayNodeOrig[map[e][n]] + i ) );
                  }
                }
              }
            }
          }
        }
      }
    }
    else if( fieldDesc.location == Location::Elem )
    {
      // Case of connectivity = Node and location = Elem
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Node;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Node
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;
            globalIndex_array const & indexArrayElem = subRegion->getReference<globalIndex_array>( fieldDesc.key );

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().nodeListString )->reference();

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  {
                    pairs.push_back( std::make_pair( indexArrayNode[map[e][n]],
                                                     fieldDesc.numComponents * indexArrayElem[e] + i ) );
                  }
                }
              }
            }
          }
        }
      }
    }
    else
    {
      // Case of connectivity = Node and location = Face
      FieldDescription fieldTmp;

      // Create a name decoration with all region names
      string nameDecoration = "";
      for( localIndex i = 0 ; i < fieldDesc.regionNames.size() ; ++i )
      {
        if( activeRegions[i] >= 0 )
        {
          nameDecoration += fieldDesc.regionNames[i];
        }
      }

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Node;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key + nameDecoration;

      createIndexArray_NodeOrFaceVersion( fieldTmp, activeRegions );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldTmp.key );

      globalIndex_array indexArrayFaceMarker( fieldDesc.numGlobalRows );
      indexArrayFaceMarker = 1;

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Face and Location::Node
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
      {
        if( activeRegions[er] >= 0 )
        {
          for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
          {
            CellElementSubRegion const * const
            subRegion = fieldDesc.regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
            integer_array const & ghostRank = subRegion->m_ghostRank;

            localIndex_array2d const &
            map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference();
            localIndex_array nodeIndices;

            for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            {
              if( ghostRank[e] < 0 )
              {
                for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                {
                  globalIndex faceId = indexArrayFace[map[e][n]];
                  subRegion->GetFaceNodes( e, n, nodeIndices );
                  if( indexArrayFaceMarker( faceId ) >= 0 )
                  {
                    // Mark this face as already visited
                    indexArrayFaceMarker( faceId ) = -1;
                    for( localIndex j = 0 ; j < nodeIndices.size() ; ++j )
                    {
                      for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                      {
                        pairs.push_back( std::make_pair( indexArrayNode[nodeIndices[j]],
                                                         fieldDesc.numComponents * faceId + i ) );
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Check if there is at least one entry
  bool emptyPairs = ( pairs.size() == 0 );

  if( !emptyPairs )
  {
    // Sort the pairs
    sort( pairs.begin(), pairs.end(), pairComparison() );

    // Remove duplicates (if someone is present)
    array1d<indexPair>::const_iterator endPairs = std::unique( pairs.begin(), pairs.end() );
    pairs.resize( endPairs - pairs.begin() );
  }

  localIndex nRowsLoc = ( emptyPairs ) ? 0 : pairs.back().first;
  localIndex nColsLoc = ( emptyPairs ) ?
                        0 :
                        std::max_element( pairs.begin(), pairs.end(), pairSecondComparison() )->second ;

  // Find the global number of rows
  localIndex_array localGather;
  CommunicationTools::allGather( nRowsLoc, localGather );
  localIndex nRows = localGather[0];
  for( localIndex p = 1 ; p < mpiSize ; ++p )
  {
    nRows = ( localGather[p] > nRows ) ? localGather[p] : nRows;
  }
  ++nRows;

  // Find the global number of columns
  CommunicationTools::allGather( nColsLoc, localGather );
  localIndex nCols = localGather[0];
  for( localIndex p = 1 ; p < mpiSize ; ++p )
  {
    nCols = ( localGather[p] > nCols ) ? localGather[p] : nCols;
  }
  ++nCols;

  // From the vector of pairs form a sparsity pattern
  if( !emptyPairs )
  {
    vectorOfPairsToCSR( pairs, nRows, nCols, connLocPatt );
  }
  else
  {
    connLocPatt.nRows = nRows;
    connLocPatt.nCols = nCols;
    connLocPatt.rowLengths.resize( nRows + 1 );
    connLocPatt.rowLengths = 0;
  }
}

// Convert a COO matrix in CSR format
void DofManager::vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                                     localIndex const nRows,
                                     localIndex const nCols,
                                     Dof_SparsityPattern & pattern ) const
{
  // Number of entries
  localIndex nnz = pairs.size();

  // Set dimensions
  pattern.nRows = nRows;
  pattern.nCols = nCols;

  // Allocate matrix with right sizes
  pattern.rowLengths.resize( nRows + 1 );
  pattern.colIndices.resize( nnz );

  // Convert
  localIndex irow0 = 0;
  localIndex irow1;
  localIndex k = 0;
  pattern.rowLengths[0] = 0;
  for( localIndex i = 0 ; i < nnz ; ++i )
  {
    irow1 = pairs[i].first;
    pattern.colIndices[k++] = pairs[i].second;
    if( irow1 > irow0 )
    {
      for( localIndex j = irow0 ; j < irow1 ; ++j )
      {
        pattern.rowLengths[j + 1] = i;
      }
      irow0 = irow1;
    }
  }
  ++irow0;
  // Add final entries to rowLengths
  for( localIndex j = irow0 ; j <= nRows ; ++j )
  {
    pattern.rowLengths[j] = nnz;
  }
}

// Perform the matrix-matrix product A*src = dst.
void DofManager::MatrixMatrixMultiply( EpetraMatrix const &A,
                                       bool const transA,
                                       EpetraMatrix const &B,
                                       bool const transB,
                                       EpetraMatrix &C,
                                       bool const call_FillComplete ) const
{
  int
  err = EpetraExt::MatrixMatrix::Multiply( *A.unwrappedPointer(), transA,
                                           *B.unwrappedPointer(), transB, *C.unwrappedPointer(), call_FillComplete );

  GEOS_ERROR_IF( err != 0, "Error thrown in matrix/matrix multiply routine" );

  // Using "call_FillComplete_on_result = false" with rectangular matrices because in this
  // case the function does not work. After the multiplication is performed, call close().
  if( !call_FillComplete )
  {
    C.close();
  }
}

// Map a global row index to local row index
localIndex DofManager::ParallelMatrixGetLocalRowID( EpetraMatrix const &A, globalIndex const index ) const
{
  return A.unwrappedPointer()->LRID( index );
}

// Return the local number of columns on each processor
// NOTE: direct use of NumMyCols() counts also for overlays. To avoid those, DomainMap() is needed
localIndex DofManager::numMyCols( EpetraMatrix const &A ) const
{
  return A.unwrappedPointer()->DomainMap().NumMyElements();
}

// Release internal memory
void DofManager::cleanUp()
{
  localIndex numFields = m_fields.size();

  // Release memory related to field descrition
  for( localIndex i = 0 ; i < numFields ; ++i )
  {
    m_fields[i].regionNames.clear();
    m_fields[i].regionPtrs.clear();
    delete m_fields[i].connLocPattern;
  }
  m_fields.clear();

  // Release memory related to connectivity description
  m_connectivity.clear();

  // Release memory related to sparsity pattern description
  for( localIndex i = 0 ; i < numFields ; ++i )
  {
    for( localIndex j = 0 ; j < numFields ; ++j )
    {
      if( m_sparsityPattern[i][j].first != nullptr )
      {
        delete m_sparsityPattern[i][j].first;
        delete m_sparsityPattern[i][j].second;
      }
    }
  }
  m_sparsityPattern.clear();
}

// Print a CSR pattern on file
void DofManager::printConnectivityLocationPattern( string const & field, string const & fileName ) const
{
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( field ),
                 "printConnectivityLocationPattern: requested field name must be already existing." );

  // get field index
  localIndex fieldIdx = fieldIndex( field );

  // Retrieve right sparsity pattern
  ParallelMatrix const * pattern = m_fields[fieldIdx].connLocPattern;

  string name;
  if( fileName.length() == 0 )
  {
    name = "pattern_" + field + ".mtx";
  }
  else
  {
    name = fileName;
  }

  // Print the selected pattern
  printParallelMatrix( *pattern, name );
}

// Print the given parallel matrix in Matrix Market format (MTX file)
void DofManager::printParallelMatrix( ParallelMatrix const & matrix, string const & fileName ) const
{
  // Ensure the ".mtx" extension
  string name( fileName );
  if( fileName.substr( fileName.find_last_of( "." ) + 1 ) != "mtx" )
  {
    name = fileName.substr( 0, fileName.find_last_of( "." ) ) + ".mtx";
  }

  EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *matrix.unwrappedPointer() );
}

// Print a CSR pattern on file or on screen
void DofManager::printSparsityPattern( Dof_SparsityPattern const & pattern, string const & fileName ) const
{
  if( fileName.length() == 0 )
  {
    // If on screen, only processor 0 write
    if( mpiRank == 0 )
    {
      for( localIndex i = 0 ; i < pattern.nRows ; ++i )
      {
        for( localIndex j = pattern.rowLengths[i] ; j < pattern.rowLengths[i + 1] ; ++j )
        {
          std::cout << i + 1 << " " << pattern.colIndices[j] + 1 << std::endl;
        }
      }
    }
  }
  else
  {
    // If on file, it is assumed that fileName take into account mpiRank
    std::ofstream fid;
    fid.open( fileName, std::ofstream::out );
    fid << "# " << pattern.nRows << " " << pattern.nCols << std::endl;
    for( localIndex i = 0 ; i < pattern.nRows ; ++i )
    {
      for( localIndex j = pattern.rowLengths[i] ; j < pattern.rowLengths[i + 1] ; ++j )
      {
        fid << i + 1 << " " << pattern.colIndices[j] + 1 << " " << 1 << std::endl;
      }
    }
    fid.close();
  }
}

// Print the coupling table on screen
void DofManager::printConnectivityMatrix() const
{
  if( mpiRank == 0 )
  {
    localIndex numFields = m_fields.size();

    std::cout << std::endl;
    for( localIndex i = 0 ; i < numFields ; ++i )
    {
      for( localIndex j = 0 ; j < numFields ; ++j )
      {
        switch( m_connectivity[i][j] )
        {
          case Connectivity::Elem:
            std::cout << " E ";
            break;
          case Connectivity::Face:
            std::cout << " F ";
            break;
          case Connectivity::Node:
            std::cout << " N ";
            break;
          case Connectivity::USER_DEFINED:
            std::cout << " U ";
            break;
          case Connectivity::None:
            std::cout << "   ";
            break;
        }
        if( j < numFields - 1 )
        {
          std::cout << "|";
        }
      }
      std::cout << std::endl;
      if( i < numFields - 1 )
      {
        for( localIndex j = 0 ; j < numFields - 1 ; ++j )
        {
          std::cout << "---|";
        }
        std::cout << "---";
      }
      std::cout << std::endl;
    }
  }
}

}
