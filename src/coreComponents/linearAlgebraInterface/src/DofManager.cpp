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

// .... COMM TOOLS :: ALL GATHER
//      TODO: move to CommunicationTools

namespace CommTools
{
void allGather( localIndex const myValue, localIndex_array & allValues )
{
#ifdef GEOSX_USE_MPI
  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
  int mpiSize;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );

  allValues.resize( mpiSize );
  array1d<int> tmpArray( mpiSize );

  int tmpValue = integer_conversion<int>( myValue );

  MPI_Allgather( &tmpValue, 1, MPI_INT, tmpArray.data(), 1, MPI_INT, MPI_COMM_GEOSX );

  for( localIndex i = 0 ; i < tmpArray.size() ; ++i )
    allValues[i] = integer_conversion<localIndex>( tmpArray[i] );
#else
  int mpiRank = 0;
  allValues.resize(1);
  allValues[0] = myValue;
#endif
}
}

// .... DOF MANAGER :: CONSTRUCTOR

DofManager::DofManager()
{
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  // we pre-allocate an oversized array to store connectivity type
  // instead of resizing it dynamically as fields are added.

  m_connectivity.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );

  for( localIndex i = 0 ; i < MAX_NUM_FIELDS ; ++i )
    for( localIndex j = 0 ; j < MAX_NUM_FIELDS ; ++j )
      m_connectivity[i][j] = Connectivity::None;

  // we pre-allocate an oversized array to store sparsity pattern type
  // instead of resizing it dynamically as fields are added.

  m_sparsityPattern.resize( MAX_NUM_FIELDS, MAX_NUM_FIELDS );
}

// .... DOF MANAGER :: SET MESH

void DofManager::setMesh( DomainPartition * const domain,
                          localIndex const meshLevelIndex,
                          localIndex const meshBodyIndex )
{
  GEOS_ERROR_IF( m_meshLevel != nullptr, "A mesh is already assigned to this DofManager." );
  m_domain = domain;
  m_meshLevel = m_domain->getMeshBodies()->GetGroup<MeshBody>( meshBodyIndex )->getMeshLevel( meshLevelIndex );
  ;
}

// .... DOF MANAGER :: FIELD INDEX

localIndex DofManager::fieldIndex( string const & key ) const
                                   {
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
    if( m_fields[i].name == key )
      return i;
  GEOS_ERROR( "Field's string key not found in list of active fields." );
  return -1;
}

// .... DOF MANAGER :: KEY IN USE

bool DofManager::keyInUse( string const & key ) const
                           {
  for( localIndex i = 0 ; i < m_fields.size() ; ++i )
    if( m_fields[i].name == key )
      return true;
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
      sumGlobalDofs += m_fields[i].numGlobalRows;
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
      sumLocalDofs += m_fields[i].numLocalRows;
    return sumLocalDofs;
  }
}

// .... DOF MANAGER :: ADD FIELD

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
  description.regionNames = regions; // TODO: if regions is empty, add list of all regions
  description.key = field + "_dof_indices";
  description.docstring = field + " dof indices";

  if( components > 1 )
    description.docstring += " (with " + std::to_string( components ) + "-component blocks)";

  m_fields.push_back( description );

  localIndex numFields = m_fields.size();
  GEOS_ERROR_IF( numFields > MAX_NUM_FIELDS, "Limit on DofManager's MAX_NUM_FIELDS exceeded." );

  // temp reference to last field

  FieldDescription & last = m_fields[numFields - 1];

  // save pointers to "active" element regions

  ElementRegionManager * const elemManager = m_meshLevel->getElemManager();

  localIndex numTotalRegions = elemManager->numRegions();
  localIndex numActiveRegions = regions.size() == 0 ? numTotalRegions : regions.size();

  last.regionPtrs.resize( numActiveRegions );
  for( localIndex er = 0 ; er < numActiveRegions ; ++er )
  {
    if( numActiveRegions < numTotalRegions )
      last.regionPtrs[er] = elemManager->GetRegion( last.regionNames[er] ); // get by name
    else
      last.regionPtrs[er] = elemManager->GetRegion( er ); // get by index

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
    last.fieldOffset = 0;

  // save field's connectivity type (self-to-self)

  m_connectivity[numFields - 1][numFields - 1] = connectivity;

  // add sparsity pattern (LC matrix)
  SparsityPattern & connLocPatt = last.connLocPattern;
  addDiagSparsityPattern( connLocPatt, numFields - 1, connectivity );

  if( connectivity != Connectivity::Elem )
  {
    //AAAAAAAAAAAAAAAA
  }

  //if (field == "pressure")
  printSparsityPattern( connLocPatt, "FILETMP" + std::to_string( mpiRank ) );

  // Tranpose the LC matrix
  SparsityPattern pattTransp;
  matrixTranpose( connLocPatt, pattTransp );

  // local pointer to the right pattern
  SparsityPattern & pattern = m_sparsityPattern[numFields - 1][numFields - 1];

  // Compute the sparsity pattern (LC'*LC)
  localIndex ierr = matrixMatrixProduct( pattTransp, connLocPatt, pattern );

  // Check for error while computing matrix matrix product
  GEOS_ERROR_IF( ierr, "addDiagSparsityPattern: error while computing the product of sparse matrices." );

  // log some basic info

  GEOS_LOG_RANK_0( "DofManager :: Added field .... " << last.docstring );
  GEOS_LOG_RANK_0( "DofManager :: Global dofs .... " << last.numGlobalRows );
  GEOS_LOG_RANK_0( "DofManager :: Field offset ... " << last.fieldOffset );
}

// .... DOF MANAGER :: CREATE INDEX ARRAY

void DofManager::createIndexArray_NodeOrFaceVersion( FieldDescription & field ) const
                                                     {
  // step 0. register an index array with default = -1

  ObjectManagerBase * baseManager =
      field.location == Location::Node ?
                                         static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() ) :
                                         static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );

  baseManager->RegisterViewWrapper<globalIndex_array>( field.key )->
                                                                  setApplyDefaultValue( -1 )->
                                                                                            setPlotLevel(
      dataRepository::PlotLevel::LEVEL_1 )->
                                          setDescription( field.docstring );

  globalIndex_array & indexArray = baseManager->getReference<globalIndex_array>( field.key );

  // step 1. loop over all active regions
  //         determine number of local rows
  //         and sequentially number objects

  field.numLocalNodes = 0;

  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
    for( localIndex esr = 0 ; esr < field.regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellBlockSubRegion * subRegion = field.regionPtrs[er]->GetSubRegion( esr );
      integer_array const & ghostRank = subRegion->m_ghostRank;

      localIndex_array2d const & map = (
          field.location == Location::Node ?
                                             subRegion->getWrapper<FixedOneToManyRelation>(
                                                 subRegion->viewKeys().nodeList )->reference() :
                                             subRegion->getWrapper<FixedOneToManyRelation>(
                                                 subRegion->viewKeys().faceList )->reference() );

      // Set which process owns the boundary nodes/faces
      for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
        if( !( ghostRank[e] < 0 ) and ghostRank[e] < mpiRank )
          for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
            indexArray[map[e][n]] = globalIndexMax;

      for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
        if( ghostRank[e] < 0 )
        {
          for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
          {
            localIndex i = map[e][n];
            if( indexArray[i] == -1 )
            {
              indexArray[i] = field.numLocalNodes;
              field.numLocalNodes++;
            }
          }
        }
    }

  // step 2. gather row counts across ranks

  localIndex_array gather;

  CommTools::allGather( field.numLocalNodes, gather );

  field.numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
    field.numGlobalRows += gather[p];

  field.firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
    field.firstLocalRow += gather[p];

  // step 3. adjust local values to reflect processor offset

  for( localIndex n = 0 ; n < indexArray.size() ; ++n )
    if( indexArray[n] != -1 )
      indexArray[n] += field.firstLocalRow;

  // step 4. synchronize across ranks

  std::map<string, string_array> fieldNames;

  if( field.location == Location::Node )
    fieldNames["node"].push_back( field.key );
  else
    fieldNames["face"].push_back( field.key );

  CommunicationTools::SynchronizeFields(
      fieldNames, m_meshLevel,
      m_domain->getReference<array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components

  field.numLocalRows = field.numLocalNodes * field.numComponents;
  field.numGlobalRows *= field.numComponents;
  field.firstLocalRow *= field.numComponents;

  // Replace globalIndexMax with -1
  for( localIndex i = 0 ; i < indexArray.size() ; ++i )
    if( indexArray[i] == globalIndexMax )
      indexArray[i] = -1;
}

// .... DOF MANAGER :: CREATE INDEX ARRAY :: ELEMENT VERSION
//      TODO: revise to look more like node version.
//            may even be able to condense to one function.

void DofManager::createIndexArray_ElemVersion( FieldDescription & field )
{
  // step 1. loop over all active regions
  //         determine number of local rows

  field.numLocalNodes = 0;
  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
  {
    field.regionPtrs[er]->forCellBlocks( [&]( CellBlockSubRegion * const subRegion )
    {
      localIndex numGhost = subRegion->GetNumberOfGhosts();
      field.numLocalNodes += subRegion->size() - numGhost;
    } );
  }

  // step 2. gather row counts across ranks

  localIndex_array gather;

  CommTools::allGather( field.numLocalNodes, gather );

  field.numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
    field.numGlobalRows += gather[p];

  field.firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
    field.firstLocalRow += gather[p];

  // step 3. loop again (sequential policy)
  //         allocate the index array
  //         set unique global indices

  globalIndex const isUnset = -1;
  globalIndex count = 0;

  for( localIndex er = 0 ; er < field.regionPtrs.size() ; ++er )
  {
    for( localIndex esr = 0 ; esr < field.regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellBlockSubRegion * subRegion = field.regionPtrs[er]->GetSubRegion( esr );

      subRegion->RegisterViewWrapper<globalIndex_array>( field.key )->
                                                                    setApplyDefaultValue( isUnset )->
                                                                                                   setPlotLevel(
          dataRepository::PlotLevel::LEVEL_1 )-> // TODO: level 1 or 2?
          setDescription( field.docstring );

      globalIndex_array & indexArray = subRegion->getReference<globalIndex_array>( field.key );
      integer_array const & ghostRank = subRegion->m_ghostRank;

      GEOS_ERROR_IF( indexArray.size() != ghostRank.size(), "Mismatch in ghost rank and index array sizes." );

      for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
        if( ghostRank[elem] < 0 )
        {
          indexArray[elem] = field.firstLocalRow + count;
          count++;
        }
    }
  }

  GEOS_ERROR_IF( count != field.numLocalNodes, "Mismatch during assignment of local row indices" );

  // step 4. synchronize across ranks

  std::map<string, string_array> fieldNames;
  fieldNames["elem"].push_back( field.key );

  CommunicationTools::SynchronizeFields(
      fieldNames, m_meshLevel,
      m_domain->getReference<array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components

  field.numLocalRows = field.numLocalNodes * field.numComponents;
  field.numGlobalRows *= field.numComponents;
  field.firstLocalRow *= field.numComponents;
}

SparsityPattern const & DofManager::getSparsityPattern( string const & rowField,
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
    rowFieldIndex = -1;

  if( colField.length() > 0 )
  {
    // check if the col field name is already added
    GEOS_ERROR_IF( !keyInUse( colField ), "getSparsityPattern: requested field name must be already existing." );

    // get col field index
    colFieldIndex = fieldIndex( colField );
  }
  else
    colFieldIndex = -1;

  if( rowFieldIndex >= 0 and colFieldIndex >= 0 )
    return m_sparsityPattern[rowFieldIndex][colFieldIndex];
  else
  {
    localIndex_array rowArray;
    localIndex_array colArray;
    if( rowFieldIndex < 0 and colFieldIndex < 0 )
    {
      rowArray.resize( m_fields.size() );
      colArray.resize( m_fields.size() );
      for( localIndex i = 0 ; i < m_fields.size() ; ++i )
      {
        rowArray[i] = i;
        colArray[i] = i;
      }
    }
    else if( rowFieldIndex < 0 )
    {
      rowArray.resize( m_fields.size() );
      colArray.resize( 1 );
      for( localIndex i = 0 ; i < m_fields.size() ; ++i )
        rowArray[i] = i;
      colArray[0] = colFieldIndex;
    }
    else
    {
      rowArray.resize( 1 );
      colArray.resize( m_fields.size() );
      rowArray[0] = rowFieldIndex;
      for( localIndex i = 0 ; i < m_fields.size() ; ++i )
        colArray[i] = i;
    }
    static SparsityPattern pattern = combineCSRMatrices( rowArray, colArray );
    return pattern;
  }
}

void DofManager::addCoupling( string const & rowField,
                              string const & colField,
                              Connectivity const connectivity,
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

  // Check if row and col regions are the same
  bool areEqual = false;
  if( rowFieldRegionNames.size() == colFieldRegionNames.size() )
    areEqual = std::equal( rowFieldRegionNames.begin(),
                           rowFieldRegionNames.end(),
                           colFieldRegionNames.begin() );

  GEOS_ERROR_IF( !areEqual, "addCoupling: regions from the two fields have to be the same." );

  // save field's connectivity type (rowField to colField)

  m_connectivity[rowFieldIndex][colFieldIndex] = connectivity;

  // local pointer to the right pattern

  SparsityPattern & pattern = m_sparsityPattern[rowFieldIndex][colFieldIndex];

  if( connectivity == Connectivity::None )
  {
    // create empty matrix
    createEmptyMatrix( rowFieldDesc.numLocalRows, colFieldDesc.numLocalRows, pattern );
  }
  else
  {
    // add sparsity pattern
    addExtraDiagSparsityPattern( pattern, rowFieldIndex, colFieldIndex, connectivity );
  }

  if( symmetric )
  {
    // Set connectivity
    m_connectivity[colFieldIndex][rowFieldIndex] = connectivity;

    // local pointer to the right pattern
    SparsityPattern & patternTransp = m_sparsityPattern[colFieldIndex][rowFieldIndex];

    // transpose the pattern
    matrixTranpose( pattern, patternTransp );
  }

}

// Get global indices for dofs connected by the connector type.
void DofManager::getIndices( globalIndex_array & indices,
                             Connectivity const connectivity,
                             localIndex const region,
                             localIndex const subregion,
                             localIndex const index,
                             string const & field ) const
                             {
  if( connectivity != Connectivity::Elem )
    getIndices( indices, connectivity, index, field );
  else
  {
    // check if the field name is already added
    GEOS_ERROR_IF( !keyInUse( field ), "getIndices: requested field name must be already existing." );

    // get field index
    localIndex const fieldIdx = fieldIndex( field );

    // check connectivity
    if( m_connectivity[fieldIdx][fieldIdx] != connectivity )
      indices.resize( 0 );
    else
    {
      // get field description
      FieldDescription const & fieldDesc = m_fields[fieldIdx];

      if( region >= fieldDesc.regionPtrs.size() )
        indices.resize( 0 );
      else
      {
        if( subregion >= fieldDesc.regionPtrs[region]->numSubRegions() )
          indices.resize( 0 );
        else
        {
          globalIndex count = 0;
          for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
          {
            for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
            {
              CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
              integer_array const & ghostRank = subRegion->m_ghostRank;

              for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
                if( ghostRank[elem] < 0 )
                {
                  if( er == region && esr == subregion && index == elem )
                  {
                    localIndex rowStart = fieldDesc.connLocPattern.rowLengths[count];
                    localIndex rowEnd = fieldDesc.connLocPattern.rowLengths[count + 1];
                    localIndex nnz = rowEnd - rowStart;
                    indices.resize( nnz );
                    std::copy( fieldDesc.connLocPattern.colIndices.data( rowStart ),
                               fieldDesc.connLocPattern.colIndices.data( rowEnd ),
                               indices.begin() );
                  }
                  count++;
                }
            }
          }
        }
      }
    }
  }
}

// Get global indices for dofs connected by the connector type.
void DofManager::getIndices( globalIndex_array & indices,
                             Connectivity const connectivity,
                             localIndex const index,
                             string const & field ) const
                             {
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( field ), "getIndices: requested field name must be already existing." );

  // get field index
  const localIndex fieldIdx = fieldIndex( field );

  // check connectivity
  if( m_connectivity[fieldIdx][fieldIdx] != connectivity )
    indices.resize( 0 );
  else
  {
    // get field description
    const FieldDescription & fieldDesc = m_fields[fieldIdx];

    if( index >= fieldDesc.numLocalRows )
      indices.resize( 0 );
    else
    {
      localIndex rowStart = fieldDesc.connLocPattern.rowLengths[index];
      localIndex rowEnd = fieldDesc.connLocPattern.rowLengths[index + 1];
      localIndex nnz = rowEnd - rowStart;
      indices.resize( nnz );
      std::copy( fieldDesc.connLocPattern.colIndices.data( rowStart ),
                 fieldDesc.connLocPattern.colIndices.data( rowEnd ),
                 indices.begin() );
    }
  }
}

// Add the specified sparsity pattern between rowField and colField to m_sparsityPattern
// collection
void DofManager::addExtraDiagSparsityPattern( SparsityPattern & pattern,
                                              localIndex const & rowFieldIndex,
                                              localIndex const & colFieldIndex,
                                              Connectivity const connectivity ) const
                                              {
  // Row field description
  FieldDescription const & rowFieldDesc = m_fields[rowFieldIndex];

  // Col field description
  FieldDescription const & colFieldDesc = m_fields[colFieldIndex];

  // Compute the CL matrices for row and col fields
  SparsityPattern rowPattern, colPattern;
  addDiagSparsityPattern( rowPattern, rowFieldIndex, connectivity );
  addDiagSparsityPattern( colPattern, colFieldIndex, connectivity );

  // Transpose the row CL matrix
  SparsityPattern rowPatternTransp;
  matrixTranpose( rowPattern, rowPatternTransp );

  // Compute the sparsity pattern (LC'*LC)
  localIndex ierr = matrixMatrixProduct( rowPatternTransp, colPattern, pattern );

  // Check for error while computing matrix matrix product
  GEOS_ERROR_IF( ierr, "addExtraDiagSparsityPattern: error while computing the product of sparse matrices." );
}

// Compute the sparsity pattern of the matrix connectivity - location for the specified
// field (diagonal entry in the m_connectivity collection)
void DofManager::addDiagSparsityPattern( SparsityPattern & connLocPatt,
                                         localIndex const & fieldIdx,
                                         Connectivity const connectivity ) const
                                         {
  // get field description
  FieldDescription const & fieldDesc = m_fields[fieldIdx];

  // array to store the matrix in COO format
  array1d<indexPair> pairs;

  if( connectivity == Connectivity::None )
  {
    // Case of no connector
    addSelfConnectedFieldSparsityPattern( pairs, fieldIdx );
  }
  else if( connectivity == Connectivity::Elem )
  {
    if( fieldDesc.location == Location::Elem )
      // Case of connector = Elem and location = Elem
      addSelfConnectedFieldSparsityPattern( pairs, fieldIdx );
    else
    {
      // Case of connector = Elem and location = Node or Face
      ObjectManagerBase * baseManager =
          fieldDesc.location == Location::Node ?
                                                 static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() ) :
                                                 static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );

      globalIndex_array & indexArray = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      globalIndex elemIndex = 0;
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;

          localIndex_array2d const & map = (
              fieldDesc.location == Location::Node ?
                                                     subRegion->getWrapper<FixedOneToManyRelation>(
                                                         subRegion->viewKeys().nodeList )->reference() :
                                                     subRegion->getWrapper<FixedOneToManyRelation>(
                                                         subRegion->viewKeys().faceList )->reference() );

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            if( ghostRank[e] < 0 )
            {
              for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  pairs.push_back( std::make_pair( elemIndex, fieldDesc.numComponents * indexArray[map[e][n]] + i ) );
              ++elemIndex;
            }
        }
    }
  }
  else if( connectivity == Connectivity::Face )
  {
    if( fieldDesc.location == Location::Face )
      // Case of connector = Face and location = Face
      addSelfConnectedFieldSparsityPattern( pairs, fieldIdx );
    else if( fieldDesc.location == Location::Elem )
    {
      // Case of connector = Face and location = Elem
      FieldDescription fieldTmp;

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Face;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Face
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;
          globalIndex_array & indexArrayElem = subRegion->getReference<globalIndex_array>( fieldDesc.key );

          localIndex_array2d const & map = subRegion->getWrapper<FixedOneToManyRelation>(
              subRegion->viewKeys().faceList )->reference();

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            if( ghostRank[e] < 0 )
              for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  if( indexArrayFace[map[e][n]] > 0 )
                    pairs.push_back(
                        std::make_pair( indexArrayFace[map[e][n]], fieldDesc.numComponents * indexArrayElem[e] + i ) );
        }
    }
    else
    {
      // Case of connector = Face and location = Node
      FieldDescription fieldTmp;

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Face;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;

          localIndex_array2d const & map = subRegion->getWrapper<FixedOneToManyRelation>(
              subRegion->viewKeys().faceList )->reference();
          localIndex_array nodeIndices;

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            if( ghostRank[e] < 0 )
            {
              for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
              {
                globalIndex faceId = indexArrayFace[map[e][n]];
                subRegion->GetFaceNodes( e, n, nodeIndices );
                if( faceId >= 0 )
                {
                  // Mark this face as already visited
                  indexArrayFace[map[e][n]] *= -1;
                  for( localIndex j = 0 ; j < nodeIndices.size() ; ++j )
                    for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                      pairs.push_back(
                          std::make_pair( faceId, fieldDesc.numComponents * indexArrayNode[nodeIndices[j]] + i ) );
                }
              }
            }
        }

      // Restore original sign in vector indexArrayFace
      for( globalIndex i = 0 ; i < indexArrayFace.size() ; ++i )
        if( indexArrayFace[i] < 0 )
          indexArrayFace[i] *= -1;
    }
  }
  else
  {
    if( fieldDesc.location == Location::Node )
      // Case of connector = Node and location = Node
      addSelfConnectedFieldSparsityPattern( pairs, fieldIdx );
    else if( fieldDesc.location == Location::Elem )
    {
      // Case of connector = Node and location = Elem
      FieldDescription fieldTmp;

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Node;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Node
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;
          globalIndex_array & indexArrayElem = subRegion->getReference<globalIndex_array>( fieldDesc.key );

          localIndex_array2d const & map = subRegion->getWrapper<FixedOneToManyRelation>(
              subRegion->viewKeys().nodeList )->reference();

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            if( ghostRank[e] < 0 )
            {
              for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
                for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                  pairs.push_back(
                      std::make_pair( indexArrayNode[map[e][n]], fieldDesc.numComponents * indexArrayElem[e] + i ) );
            }
        }
    }
    else
    {
      // Case of connector = Node and location = Face
      FieldDescription fieldTmp;

      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.regionPtrs = fieldDesc.regionPtrs;
      fieldTmp.location = Location::Node;
      fieldTmp.numComponents = 1;
      fieldTmp.regionNames = fieldDesc.regionNames;
      fieldTmp.key = fieldDesc.key;

      createIndexArray_NodeOrFaceVersion( fieldTmp );

      ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getFaceManager() );
      globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      baseManager = static_cast<ObjectManagerBase*>( m_meshLevel->getNodeManager() );
      globalIndex_array & indexArrayNode = baseManager->getReference<globalIndex_array>( fieldDesc.key );

      // Compute the transpose of the sparsity pattern, i.e., Connectivity::Face and Location::Node
      for( localIndex er = 0 ; er < fieldDesc.regionPtrs.size() ; ++er )
        for( localIndex esr = 0 ; esr < fieldDesc.regionPtrs[er]->numSubRegions() ; esr++ )
        {
          CellBlockSubRegion * subRegion = fieldDesc.regionPtrs[er]->GetSubRegion( esr );
          integer_array const & ghostRank = subRegion->m_ghostRank;

          localIndex_array2d const & map = subRegion->getWrapper<FixedOneToManyRelation>(
              subRegion->viewKeys().faceList )->reference();
          localIndex_array nodeIndices;

          for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
            if( ghostRank[e] < 0 )
              for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
              {
                globalIndex faceId = indexArrayFace[map[e][n]];
                subRegion->GetFaceNodes( e, n, nodeIndices );
                if( faceId >= 0 )
                {
                  // Mark this face as already visited
                  indexArrayFace[map[e][n]] *= -1;
                  for( localIndex j = 0 ; j < nodeIndices.size() ; ++j )
                    for( localIndex i = 0 ; i < fieldDesc.numComponents ; ++i )
                      pairs.push_back(
                          std::make_pair( indexArrayNode[nodeIndices[j]], fieldDesc.numComponents * faceId + i ) );
                }
              }
        }

      // Restore original sign in vector indexArrayFace
      for( globalIndex i = 0 ; i < indexArrayFace.size() ; ++i )
        if( indexArrayFace[i] < 0 )
          indexArrayFace[i] *= -1;
    }
  }

  // Sort the pairs
  sort( pairs.begin(), pairs.end(), pairComparison() );

  // From the vector of pairs form a sparsity pattern
  localIndex nCols = integer_conversion<localIndex>( fieldDesc.numGlobalRows );
  vectorOfPairsToCSR( pairs, nCols, connLocPatt );
}

// Specialized function for cases when the connectors are the locations, i.e. there is no
// communications with other entities
void DofManager::addSelfConnectedFieldSparsityPattern( array1d<indexPair> & pairs, localIndex const & fieldIdx ) const
                                                       {
  // get field description
  FieldDescription const & fieldDesc = m_fields[fieldIdx];

  // get the local number of rows
  localIndex const localNodes = fieldDesc.numLocalNodes;

  for( localIndex i = 0 ; i < localNodes ; ++i )
    for( localIndex j = 0 ; j < fieldDesc.numComponents ; ++j )
      pairs.push_back( std::make_pair( i, fieldDesc.firstLocalRow + fieldDesc.numComponents * i + j ) );
}

// Convert a COO matrix in CSR format
void DofManager::vectorOfPairsToCSR( array1d<indexPair> const pairs, localIndex const nCols,
                                     SparsityPattern & pattern ) const
                                     {
  // Number of entries
  localIndex nnz = pairs.size();

  // Allocate matrix with right sizes
  pattern.colIndices.resize( nnz );

  // Convert
  localIndex irow0 = 0;
  localIndex irow1;
  localIndex k = 0;
  pattern.rowLengths.push_back( 0 );
  for( localIndex i = 0 ; i < nnz ; ++i )
  {
    irow1 = pairs[i].first;
    pattern.colIndices[k++] = pairs[i].second;
    if( irow1 > irow0 )
    {
      for( localIndex j = irow0 ; j < irow1 ; ++j )
        pattern.rowLengths.push_back( i );
      irow0 = irow1;
    }
  }
  // Add last entry to rowLengths
  pattern.rowLengths.push_back( nnz );

  // Set dimensions
  pattern.nRows = pattern.rowLengths.size() - 1;
  pattern.nCols = nCols;
}

// Create an empty matrix
void DofManager::createEmptyMatrix( localIndex const nRows,
                                    localIndex const nCols,
                                    SparsityPattern & pattern ) const
                                    {
  pattern.nRows = nRows;
  pattern.nCols = nCols;
  pattern.rowLengths.resize( nRows + 1 );
  pattern.rowLengths = 0;
  pattern.colIndices.resize( 0 );
}

// Combine several CSR matrices according to the m_sparsityPattern 2d structure
SparsityPattern DofManager::combineCSRMatrices( localIndex_array const rowArray,
                                                localIndex_array const colArray ) const
                                                {
  // Compute number of rows and columns and shifts (of rows and columns) for the global matrix
  localIndex nRowsTot = 0;
  localIndex nColsTot = 0;
  localIndex nnzTot = 0;
  localIndex_array rowShift( rowArray.size() );
  rowShift[0] = 0;
  localIndex_array colShift( colArray.size() );
  colShift[0] = 0;
  for( localIndex iGlo = 0 ; iGlo < rowArray.size() ; ++iGlo )
  {
    // Retrieve size from the diagonal block
    nRowsTot += m_fields[iGlo].numGlobalRows;
    for( localIndex jGlo = 0 ; jGlo < colArray.size() ; ++jGlo )
    {
      nnzTot += m_sparsityPattern[iGlo][jGlo].colIndices.size();
      if( iGlo == 0 and jGlo < colArray.size() - 1 )
        colShift[jGlo + 1] = colShift[jGlo] + m_fields[jGlo].numGlobalRows;
      if( iGlo == 0 )
        // Retrieve size from the diagonal block
        nColsTot += m_fields[jGlo].numGlobalRows;
    }
    if( iGlo < rowArray.size() - 1 )
      rowShift[iGlo + 1] = rowShift[iGlo] + m_fields[iGlo].numGlobalRows;
  }

  // Allocate global matrix
  SparsityPattern pattern;
  pattern.rowLengths.resize( nRowsTot + 1 );
  pattern.colIndices.resize( nnzTot );
  pattern.nRows = nRowsTot;
  pattern.nCols = nColsTot;

  localIndex_array & rowPtr = pattern.rowLengths;
  globalIndex_array & colIdx = pattern.colIndices;
  rowPtr[0] = 0;

  // Copy local matrices in the right position
  globalIndex k = 0;
  for( localIndex iGlo = 0 ; iGlo < rowArray.size() ; ++iGlo )
  {
    localIndex nLocRows = m_fields[iGlo].numGlobalRows;
    for( localIndex iLoc = 0 ; iLoc < nLocRows ; ++iLoc )
    {
      for( localIndex jGlo = 0 ; jGlo < colArray.size() ; ++jGlo )
      {
        SparsityPattern const & patternLoc = m_sparsityPattern[iGlo][jGlo];
        if( patternLoc.nRows > 0 )
        {
          localIndex_array const & rowPtrLoc = patternLoc.rowLengths;
          globalIndex_array const & colIdxLoc = patternLoc.colIndices;
          for( localIndex j = rowPtrLoc[iLoc] ; j < rowPtrLoc[iLoc + 1] ; ++j )
            colIdx[k++] = colIdxLoc[j] + colShift[jGlo];
        }
      }
      rowPtr[rowShift[iGlo] + iLoc + 1] = k;
    }
  }

  return pattern;
}

// Matrix-matrix product
localIndex DofManager::matrixMatrixProduct( SparsityPattern const & matA,
                                            SparsityPattern const & matB,
                                            SparsityPattern & matC ) const
                                            {
  localIndex ierr = 0;

  // check for error
  if( matA.nCols != matB.nRows )
  {
    ierr = 1;
    return ierr;
  }

  // Retrieve pointers
  localIndex_array const & rowPtrA = matA.rowLengths;
  globalIndex_array const & colIdxA = matA.colIndices;
  localIndex_array const & rowPtrB = matB.rowLengths;
  globalIndex_array const & colIdxB = matB.colIndices;
  localIndex_array & rowPtrC = matC.rowLengths;
  globalIndex_array & colIdxC = matC.colIndices;

  // Allocate scratch arrays:
  // - JWN1 collects the row of the product;
  // - WN1 marks if a column index is already present
  globalIndex_array JWN1( matB.nCols );
  localIndex_array WN1( matB.nCols );
  // Initialize scratch array
  WN1 = -1;

  // Set nRows and nCols for C
  matC.nRows = matA.nRows;
  matC.nCols = matB.nCols;

  // Allocate rowLengths for the product matrix
  rowPtrC.resize( matC.nRows + 1 );

  // Estimated size for colIdxC
  globalIndex nnzC = 20 * matC.nRows;
  colIdxC.resize( nnzC );

  rowPtrC[0] = 0;
  for( localIndex i = 0 ; i < matA.nRows ; ++i )
  {
    globalIndex ind = 0;
    for( localIndex j = rowPtrA[i] ; j < rowPtrA[i + 1] ; ++j )
    {
      globalIndex idColA = colIdxA[j];
      localIndex startRowB = rowPtrB[idColA];
      localIndex endRowB = rowPtrB[idColA + 1];
      // Combination of columns from matrix A
      for( localIndex k = startRowB ; k < endRowB ; ++k )
      {
        globalIndex ii = colIdxB[k];
        if( WN1[ii] == -1 )
        {
          JWN1[ind] = ii;
          WN1[ii] = ind;
          ind++;
        }
      }
    }
    rowPtrC[i + 1] = rowPtrC[i] + ind;
    std::sort( JWN1.begin(), JWN1.data( ind ) );
    if( rowPtrC[i + 1] > nnzC )
    {
      nnzC = rowPtrC[i + 1] > 2 * nnzC ? rowPtrC[i + 1] : 2 * nnzC;
      colIdxC.resize( nnzC );
    }
    std::copy( JWN1.begin(), JWN1.data( ind ), colIdxC.data( rowPtrC[i] ) );
    for( globalIndex j = 0 ; j < ind ; ++j )
      WN1[JWN1[j]] = -1;
  }
  colIdxC.resize( rowPtrC[matC.nRows] );

  return ierr;
}

void DofManager::matrixTranpose( SparsityPattern const & matA,
                                 SparsityPattern & matB ) const
                                 {
  // Retrieve pointers
  localIndex_array const & rowPtrA = matA.rowLengths;
  globalIndex_array const & colIdxA = matA.colIndices;
  localIndex_array & rowPtrB = matB.rowLengths;
  globalIndex_array & colIdxB = matB.colIndices;

  // Set nRows and nCols for B
  matB.nRows = matA.nCols;
  matB.nCols = matA.nRows;

  // Allocate memory for the matrix transpose
  localIndex nnz = colIdxA.size();
  rowPtrB.resize( matB.nRows + 1 );
  colIdxB.resize( nnz );

  // Allocate scratch to count nnz for B rows
  localIndex_array WN1( matB.nRows );

  // Explore A by rows counting the nnz of each column
  rowPtrB = 0;
  for( localIndex i = 0 ; i < matA.nRows ; ++i )
    for( localIndex j = rowPtrA[i] ; j < rowPtrA[i + 1] ; ++j )
      rowPtrB[colIdxA[j]] += 1;

  // Create rowLengths for matrix B
  WN1[0] = 0;
  for( localIndex i = 1 ; i < matA.nCols ; ++i )
    WN1[i] = WN1[i - 1] + rowPtrB[i - 1];

  // Copy WN1 into rowPtrB
  std::copy( WN1.begin(), WN1.end(), rowPtrB.begin() );

  // Set the last entry (nnz)
  rowPtrB[matA.nCols] = nnz;

  // Transpose A
  for( localIndex i = 0 ; i < matA.nRows ; ++i )
    for( localIndex j = rowPtrA[i] ; j < rowPtrA[i + 1] ; ++j )
      colIdxB[WN1[colIdxA[j]]++] = i;
}

// Print a CSR pattern on file
void DofManager::printSparsityPattern( string const & field, string fileName ) const
                                       {
  // check if the field name is already added
  GEOS_ERROR_IF( !keyInUse( field ), "printSparsityPattern: requested field name must be already existing." );

  // get field index
  localIndex fieldIdx = fieldIndex( field );

  // Retrieve right sparsity pattern
  SparsityPattern & pattern = m_sparsityPattern[fieldIdx][fieldIdx];

  if( fileName.length() == 0 )
    fileName = "pattern_" + field + ".csr";

  // Print the selected pattern
  printSparsityPattern( pattern, fileName );
}

// Print a CSR pattern on file
void DofManager::printCoupledSparsityPattern( string const & rowField, string const & colField, string fileName ) const
                                              {
  // check if the row field name is already added
  GEOS_ERROR_IF( !keyInUse( rowField ), "printSparsityPattern: requested field name must be already existing." );

  // check if the col field name is already added
  GEOS_ERROR_IF( !keyInUse( colField ), "printSparsityPattern: requested field name must be already existing." );

  // get row field index
  localIndex rowFieldIndex = fieldIndex( rowField );

  // get column field index
  localIndex colFieldIndex = fieldIndex( colField );

  // Retrieve right sparsity pattern
  SparsityPattern & pattern = m_sparsityPattern[rowFieldIndex][colFieldIndex];

  if( fileName.length() == 0 )
    fileName = "pattern_" + rowField + "_" + colField + "_" + std::to_string( mpiRank ) + ".csr";

  // Print the selected pattern
  printSparsityPattern( pattern, fileName );
}

// Print a CSR pattern on file or on screen
void DofManager::printSparsityPattern( SparsityPattern const & pattern, string const & fileName ) const
                                       {
  if( fileName.length() == 0 )
  {
    // If on screen, only processor 0 write
    if( mpiRank == 0 )
      for( localIndex i = 0 ; i < pattern.nRows ; ++i )
        for( localIndex j = pattern.rowLengths[i] ; j < pattern.rowLengths[i + 1] ; ++j )
          std::cout << i + 1 << " " << pattern.colIndices[j] + 1 << std::endl;
  }
  else
  {
    // If on file, it is assumed that fileName take into account mpiRank
    std::ofstream fid;
    fid.open( fileName, std::ofstream::out );
    fid << "# " << pattern.nRows << " " << pattern.nCols << std::endl;
    for( localIndex i = 0 ; i < pattern.nRows ; ++i )
      for( localIndex j = pattern.rowLengths[i] ; j < pattern.rowLengths[i + 1] ; ++j )
        fid << i + 1 << " " << pattern.colIndices[j] + 1 << " " << 1 << std::endl;
    fid.close();
  }
}

// Print the coupling table on screen
void DofManager::printCoupling() const
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
          case Connectivity::None:
            std::cout << "   ";
            break;
        }
        if( j < numFields - 1 )
          std::cout << "|";
      }
      std::cout << std::endl;
      if( i < numFields - 1 )
      {
        for( localIndex j = 0 ; j < numFields - 1 ; ++j )
          std::cout << "---|";
        std::cout << "---";
      }
      std::cout << std::endl;
    }
  }
}

}
