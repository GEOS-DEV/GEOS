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
 * @file createConnLocPattern.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_TESTS_CREATECONNLOCPATTERN_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_TESTS_CREATECONNLOCPATTERN_HPP_

//#include "physicsSolvers/SimpleSolvers/LaplaceFEM.hpp"
#include "DofManager.hpp"

namespace geosx
{

/**
 * Definifion for entries of sparse matrix in COO format
 */
typedef std::pair<localIndex, globalIndex> indexPair;

/**
 * Compare structure used to create CSR matrix from COO format
 */
struct pairComparison
{
  inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
  {
    if( lhs.first < rhs.first )
      return true;
    else if( lhs.first == rhs.first )
      return lhs.second < rhs.second;
    else
      return false;
  }
};

/**
 * Compare second element of a pair
 */
struct pairSecondComparison
{
  inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
  {
    return ( lhs.second < rhs.second );
  }
};

/**
 * Create face index array
 */
void createIndexArray_FaceVersion( DomainPartition * const domain,
                                   MeshLevel * const meshLevel,
                                   const string & fieldName,
                                   array1d<ElementRegion*> const & regionPtrs,
                                   localIndex const numComponents,
                                   localIndex & numLocalNodes,
                                   localIndex & numLocalRows,
                                   globalIndex & numGlobalRows,
                                   globalIndex & firstLocalRow,
                                   globalIndex & firstLocalConnectivity );

/**
 * Create element index array
 */
void createIndexArray_ElemVersion( DomainPartition * const domain,
                                   MeshLevel * const meshLevel,
                                   const string & fieldName,
                                   array1d<ElementRegion*> const & regionPtrs,
                                   localIndex const numComponents,
                                   localIndex & numLocalNodes,
                                   localIndex & numLocalRows,
                                   globalIndex & numGlobalRows,
                                   globalIndex & firstLocalRow,
                                   globalIndex & firstLocalConnectivity );

/**
 * Convert a COO matrix in CSR format
 */
void vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                         localIndex const nRows,
                         localIndex const nCols,
                         Dof_SparsityPattern & pattern );

/**
 * Case of connectivity = Face and location = Elem
 */
void addDiagSparsityPattern( Dof_SparsityPattern & connLocPatt,
                             DomainPartition * domain,
                             MeshLevel * meshLevel,
                             const string & fieldName,
                             array1d<ElementRegion*> const & regionPtrs,
                             localIndex const numComponents );

/**
 * Create the connectivity-location pattern for the TPFA finite volume approach
 */
void createConnLocPattern( DomainPartition * const domain,
                           localIndex const meshBodyIndex,
                           localIndex const meshLevelIndex,
                           localIndex const numComponents,
                           ParallelMatrix & connLocPattern );

///**
// * Create the location-location pattern for the Laplace equation with FEM
// */
//void createLocLocPattern( DomainPartition * const domain,
//                          localIndex const meshBodyIndex,
//                          localIndex const meshLevelIndex,
//                          ParallelMatrix & locLocPattern );

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_TESTS_CREATECONNLOCPATTERN_HPP_ */
