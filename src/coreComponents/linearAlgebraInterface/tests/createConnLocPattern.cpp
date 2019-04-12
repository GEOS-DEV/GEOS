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
 * @file createConnLocPattern.cpp
 */

#include "createConnLocPattern.hpp"

namespace geosx
{

// Create face index array
void createIndexArray_FaceVersion( DomainPartition * const domain,
                                   MeshLevel * const meshLevel,
                                   const string & fieldName,
                                   array1d<ElementRegion*> const & regionPtrs,
                                   localIndex const numComponents,
                                   localIndex & numLocalNodes,
                                   localIndex & numLocalRows,
                                   globalIndex & numGlobalRows,
                                   globalIndex & firstLocalRow,
                                   globalIndex & firstLocalConnectivity )
{
  int mpiSize;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  // step 0. register an index array with default = LocationStatus::notAssigned
  ObjectManagerBase *
  baseManager = static_cast<ObjectManagerBase*>( meshLevel->getFaceManager() );

  baseManager->RegisterViewWrapper<globalIndex_array>( fieldName )->
    setApplyDefaultValue( static_cast<globalIndex>( DofManager::LocationStatus::notAssigned ) )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );

  globalIndex_array & indexArray = baseManager->getReference<globalIndex_array>( fieldName );

  // step 1. loop over all active regions
  //         determine number of local rows
  //         and sequentially number objects
  numLocalNodes = 0;
  localIndex numLocalNodesWithGhost = 0;
  globalIndex numLocalConnectivity = 0;

  for( localIndex er = 0 ; er < regionPtrs.size() ; ++er )
  {
    for( localIndex esr = 0 ; esr < regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellElementSubRegion const * const subRegion = regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
      integer_array const & ghostRank = subRegion->m_ghostRank;

      localIndex_array2d const &
      map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference() ;

      // Set which process owns the boundary nodes/faces
      for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
      {
        if( !( ghostRank[e] < 0 ) and ghostRank[e] < mpiRank )
        {
          for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
          {
            indexArray[map[e][n]] = static_cast<globalIndex>( DofManager::LocationStatus::notMyGhostLocation );
          }
        }
      }

      for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
      {
        for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
        {
          localIndex i = map[e][n];
          if( indexArray[i] == static_cast<globalIndex>( DofManager::LocationStatus::notAssigned ) )
          {
            indexArray[i] = numLocalNodesWithGhost;
            numLocalNodesWithGhost++;
            if( ghostRank[e] < 0 )
            {
              ++numLocalNodes;
            }
          }
        }
        ++numLocalConnectivity;
      }
    }
  }

  // step 2. gather row counts across ranks
  localIndex_array localGather;

  CommunicationTools::allGather( numLocalNodes, localGather );

  numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    numGlobalRows += localGather[p];
  }

  firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    firstLocalRow += localGather[p];
  }

  // for starting the connectivity numbering
  globalIndex_array globalGather;
  CommunicationTools::allGather( numLocalConnectivity, globalGather );

  firstLocalConnectivity = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    firstLocalConnectivity += globalGather[p];
  }

  // step 3. adjust local values to reflect processor offset
  for( localIndex n = 0 ; n < indexArray.size() ; ++n )
  {
    if( indexArray[n] != static_cast<globalIndex>( DofManager::LocationStatus::notAssigned ) )
    {
      indexArray[n] += firstLocalRow;
    }
  }

  // step 4. synchronize across ranks
  std::map<string, string_array> fieldNames;
  fieldNames["face"].push_back( fieldName );

  CommunicationTools::
  SynchronizeFields( fieldNames, meshLevel,
                     domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components
  numLocalRows = numLocalNodes * numComponents;
  numGlobalRows *= numComponents;
  firstLocalRow *= numComponents;

  // Replace LocationStatus::notMyGhostLocation with LocationStatus::notAssigned (if any)
  for( localIndex i = 0 ; i < indexArray.size() ; ++i )
  {
    if( indexArray[i] == static_cast<globalIndex>( DofManager::LocationStatus::notMyGhostLocation ) )
    {
      indexArray[i] = static_cast<globalIndex>( DofManager::LocationStatus::notAssigned );
    }
  }
}

// Create element index array
void createIndexArray_ElemVersion( DomainPartition * const domain,
                                   MeshLevel * const meshLevel,
                                   const string & fieldName,
                                   array1d<ElementRegion*> const & regionPtrs,
                                   localIndex const numComponents,
                                   localIndex & numLocalNodes,
                                   localIndex & numLocalRows,
                                   globalIndex & numGlobalRows,
                                   globalIndex & firstLocalRow,
                                   globalIndex & firstLocalConnectivity )
{
  int mpiSize;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  // step 1. loop over all active regions
  //         determine number of local rows
  numLocalNodes = 0;
  for( localIndex er = 0 ; er < regionPtrs.size() ; ++er )
  {
    regionPtrs[er]->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion const * const subRegion )
    {
      localIndex numGhost = subRegion->GetNumberOfGhosts();
      numLocalNodes += subRegion->size() - numGhost;
    } );
  }

  // step 2. gather row counts across ranks
  localIndex_array localGather;

  CommunicationTools::allGather( numLocalNodes, localGather );

  numGlobalRows = 0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    numGlobalRows += localGather[p];
  }

  firstLocalRow = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    firstLocalRow += localGather[p];
  }

  // step 3. loop again (sequential policy)
  //         allocate the index array
  //         set unique global indices
  globalIndex count = 0;

  for( localIndex er = 0 ; er < regionPtrs.size() ; ++er )
  {
    for( localIndex esr = 0 ; esr < regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellElementSubRegion * const subRegion = regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );

      subRegion->RegisterViewWrapper<globalIndex_array>( fieldName )->
        setApplyDefaultValue( static_cast<globalIndex>( DofManager::LocationStatus::notAssigned ) )->
        setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );

      globalIndex_array & indexArray = subRegion->getReference<globalIndex_array>( fieldName );
      integer_array const & ghostRank = subRegion->m_ghostRank;

      GEOS_ERROR_IF( indexArray.size() != ghostRank.size(), "Mismatch in ghost rank and index array sizes." );

      for( localIndex elem = 0 ; elem < ghostRank.size() ; ++elem )
      {
        if( ghostRank[elem] < 0 )
        {
          indexArray[elem] = firstLocalRow + count;
          count++;
        }
      }
    }
  }

  GEOS_ERROR_IF( count != numLocalNodes, "Mismatch during assignment of local row indices" );

  // for starting the connectivity numbering
  globalIndex_array globalGather;

  CommunicationTools::allGather( count, globalGather );

  firstLocalConnectivity = 0;
  for( localIndex p = 0 ; p < mpiRank ; ++p )
  {
    firstLocalConnectivity += globalGather[p];
  }

  // step 4. synchronize across ranks
  std::map<string, string_array> fieldNames;
  fieldNames["elem"].push_back( fieldName );

  CommunicationTools::
  SynchronizeFields( fieldNames, meshLevel,
                     domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components
  numLocalRows = numLocalNodes * numComponents;
  numGlobalRows *= numComponents;
  firstLocalRow *= numComponents;
}

// Convert a COO matrix in CSR format
void vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                         localIndex const nRows,
                         localIndex const nCols,
                         Dof_SparsityPattern & pattern )
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

// Case of connectivity = Face and location = Elem
void addDiagSparsityPattern( Dof_SparsityPattern & connLocPatt,
                             DomainPartition * domain,
                             MeshLevel * meshLevel,
                             const string & fieldName,
                             array1d<ElementRegion*> const & regionPtrs,
                             localIndex const numComponents)
{

  localIndex numLocalNodes;
  localIndex numLocalRows;
  globalIndex numGlobalRows;
  globalIndex firstLocalRow;
  globalIndex firstLocalConnectivity;

  createIndexArray_FaceVersion( domain,
                                meshLevel,
                                "userDefined",
                                regionPtrs,
                                1,
                                numLocalNodes,
                                numLocalRows,
                                numGlobalRows,
                                firstLocalRow,
                                firstLocalConnectivity );

  ObjectManagerBase * baseManager = static_cast<ObjectManagerBase*>( meshLevel->getFaceManager() );
  globalIndex_array & indexArrayFace = baseManager->getReference<globalIndex_array>( fieldName );

  array1d<indexPair> pairs;

  // Compute the transpose of the sparsity pattern, i.e., Connectivity::Elem and Location::Face
  for( localIndex er = 0 ; er < regionPtrs.size() ; ++er )
  {
    for( localIndex esr = 0 ; esr < regionPtrs[er]->numSubRegions() ; esr++ )
    {
      CellElementSubRegion const * const
      subRegion = regionPtrs[er]->GetSubRegion<CellElementSubRegion>( esr );
      integer_array const & ghostRank = subRegion->m_ghostRank;
      globalIndex_array const & indexArrayElem = subRegion->getReference<globalIndex_array>( fieldName );

      localIndex_array2d const &
      map = subRegion->getWrapper<FixedOneToManyRelation>( subRegion->viewKeys().faceListString )->reference();

      for( localIndex e = 0 ; e < map.size( 0 ) ; ++e )
      {
        if( ghostRank[e] < 0 )
        {
          for( localIndex n = 0 ; n < map.size( 1 ) ; ++n )
          {
            for( localIndex i = 0 ; i < numComponents ; ++i )
            {
              if( indexArrayFace[map[e][n]] >= 0 )
              {
                pairs.push_back( std::make_pair( indexArrayFace[map[e][n]],
                                                 numComponents * indexArrayElem[e] + i ) );
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

  int mpiSize;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );

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

// Create the connectivity-location pattern for the TPFA finite volume approach
void createConnLocPattern( DomainPartition * const domain,
                           localIndex const meshBodyIndex,
                           localIndex const meshLevelIndex,
                           localIndex const numComponents,
                           ParallelMatrix & connLocPattern )
{
  MeshLevel * const
  meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>( meshBodyIndex )->
              getMeshLevel( integer_conversion<int>( meshLevelIndex ) );

  // save pointers to "active" element regions
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  auto const & regionListPtr = elemManager->GetRegions().keys();
  string_array regionNames( regionListPtr.size() );

  for( auto& regionPtr : regionListPtr )
  {
    regionNames[regionPtr.second] = regionPtr.first;
  }

  localIndex numTotalRegions = elemManager->numRegions();
  localIndex numActiveRegions = regionNames.size();

  array1d<ElementRegion*> regionPtrs( numActiveRegions );
  for( localIndex er = 0 ; er < numActiveRegions ; ++er )
  {
    // Get region by name
    regionPtrs[er] = elemManager->GetRegion( regionNames[er] );
    GEOS_ERROR_IF( regionPtrs[er] == nullptr, "Specified element region not found" );
  }

  localIndex numLocalNodes, numLocalRows;
  globalIndex numGlobalRows, firstLocalRow, firstLocalConnectivity;
  createIndexArray_ElemVersion( domain,
                                meshLevel,
                                "userDefined",
                                regionPtrs,
                                1,
                                numLocalNodes,
                                numLocalRows,
                                numGlobalRows,
                                firstLocalRow,
                                firstLocalConnectivity );

  // add sparsity pattern (LC matrix)
  Dof_SparsityPattern connLocPattLocal;
  addDiagSparsityPattern( connLocPattLocal,
                          domain,
                          meshLevel,
                          "userDefined",
                          regionPtrs,
                          numComponents);

  // TRILINOS interface
  localIndex maxEntriesPerRow = 0;
  for( localIndex i = 0 ; i < connLocPattLocal.nRows ; ++i )
  {
    maxEntriesPerRow = std::max( maxEntriesPerRow,
                                 connLocPattLocal.rowLengths[i + 1] - connLocPattLocal.rowLengths[i] );
  }

  connLocPattern.createWithGlobalSize( connLocPattLocal.nRows,
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
      connLocPattern.insert( i,
                             connLocPattLocal.colIndices.data( connLocPattLocal.rowLengths[i] ),
                             values.data(),
                             nnz );
    }
  }
  connLocPattern.close();

  regionPtrs.clear();
}

// Create the location-location pattern for the Laplace equation with FEM
void createLocLocPattern( DomainPartition * const domain,
                          localIndex const meshBodyIndex,
                          localIndex const meshLevelIndex,
                          ParallelMatrix & locLocPattern )
{
  LaplaceFEM laplaceFEM( "laplaceFEM", domain );

  MeshLevel * const
  mesh = domain->getMeshBodies()->GetGroup<MeshBody>( meshBodyIndex )->
         getMeshLevel( integer_conversion<int>( meshLevelIndex ) );
  NodeManager * const nodeManager = mesh->getNodeManager();

  nodeManager->RegisterViewWrapper<array1d<globalIndex> >( laplaceFEM.laplaceFEMViewKeys.blockLocalDofNumberString )->
    setApplyDefaultValue( -1 )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )->
    setDescription("Global DOF numbers for the primary field variable");

  localIndex n_ghost_rows = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows = nodeManager->size()-n_ghost_rows;
  globalIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  laplaceFEM.SetNumRowsAndTrilinosIndices( nodeManager,
                                           n_local_rows,
                                           n_global_rows,
                                           displacementIndices,
                                           0 );
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back( laplaceFEM.laplaceFEMViewKeys.blockLocalDofNumberString );

  CommunicationTools::
  SynchronizeFields(fieldNames, mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  Epetra_Map rowMap( n_global_rows, 0, Epetra_MpiComm( MPI_COMM_GEOSX ) );
  Epetra_FECrsGraph sparsity( Copy, rowMap, 0 );
  laplaceFEM.SetSparsityPattern( domain, &sparsity );
  sparsity.FillComplete();

  locLocPattern.create( sparsity );
  locLocPattern.unwrappedPointer()->PutScalar( 1.0 );
}

}
