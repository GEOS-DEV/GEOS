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

#include "DofManager.hpp"

namespace geosx
{

/**
 * @typedef indexPair
 *
 * @brief Definifion for entries of sparse matrix in COO format.
 */
typedef std::pair<localIndex, globalIndex> indexPair;

/**
 * @struct pairComparison
 *
 * @brief Compare structure used to create CSR matrix from COO format.
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
 * @struct pairSecondComparison
 *
 * @brief Compare second element of a pair.
 */
struct pairSecondComparison
{
  inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
  {
    return ( lhs.second < rhs.second );
  }
};

/**
 * @brief Create face index array.
 *
 * @param [in]  domain DomainPartition the input domain.
 * @param [in]  meshLevel MeshLevel the mesh level.
 * @param [in]  fieldName string the field name.
 * @param [in]  regionPtrs array1d<ElementRegion> array of regions where the field is defined.
 * @param [in]  numComponents localIndex number of components (for vector fields).
 * @param [out] numLocalNodes localIndex number of local nodes.
 * @param [out] numLocalRows localIndex number of local rows.
 * @param [out] numGlobalRows globalIndex number of global rows.
 * @param [out] firstLocalRow globalIndex first row on this processor (without field offset).
 * @param [out] firstLocalConnectivity globalIndex first connector on this processor.
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
 * @brief Create element index array.
 *
 * @param [in]  domain DomainPartition the input domain.
 * @param [in]  meshLevel MeshLevel the mesh level.
 * @param [in]  fieldName string the field name.
 * @param [in]  regionPtrs array1d<ElementRegion> array of regions where the field is defined.
 * @param [in]  numComponents localIndex number of components (for vector fields).
 * @param [out] numLocalNodes localIndex number of local nodes.
 * @param [out] numLocalRows localIndex number of local rows.
 * @param [out] numGlobalRows globalIndex number of global rows.
 * @param [out] firstLocalRow globalIndex first row on this processor (without field offset).
 * @param [out] firstLocalConnectivity globalIndex first connector on this processor.
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
 * @brief Convert a COO matrix in CSR format.
 *
 * @param [in]  pairs array1d<indexPair> array of pairs (irow, jcol): COO format.
 * @param [in]  nRows localIndex number of rows.
 * @param [in]  nCols localIndex number of columns.
 * @param [out] pattern Dof_SparsityPattern matrix in CSR format.
 */
void vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                         localIndex const nRows,
                         localIndex const nCols,
                         Dof_SparsityPattern & pattern );

/**
 * @brief Case of connectivity = Face and location = Elem.
 *
 * @param [out] connLocPatt Dof_SparsityPattern the local sparsity pattern.
 * @param [in]  domain DomainPartition the input domain.
 * @param [in]  meshLevel MeshLevel the mesh level.
 * @param [in]  fieldName string the field name.
 * @param [in]  regionPtrs array1d<ElementRegion> array of regions where the field is defined.
 * @param [in]  numComponents localIndex number of components (for vector fields).
 */
void addDiagSparsityPattern( Dof_SparsityPattern & connLocPatt,
                             DomainPartition * const domain,
                             MeshLevel * const meshLevel,
                             const string & fieldName,
                             array1d<ElementRegion*> const & regionPtrs,
                             localIndex const numComponents );

/**
 * @brief Create the connectivity-location pattern for the TPFA finite volume approach.
 *
 * @param [in]  domain DomainPartition the input domain.
 * @param [in]  meshLevel MeshLevel the mesh level.
 * @param [in]  fieldName string the field name.
 * @param [in]  regionPtrs array1d<ElementRegion> array of regions where the field is defined.
 * @param [in]  numComponents localIndex number of components (for vector fields).
 * @param [out] connLocPattern ParallelMatrix the output sparsity pattern.
 */
void createConnLocPattern( DomainPartition * const domain,
                           localIndex const meshBodyIndex,
                           localIndex const meshLevelIndex,
                           localIndex const numComponents,
                           ParallelMatrix & connLocPattern );

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_TESTS_CREATECONNLOCPATTERN_HPP_ */
