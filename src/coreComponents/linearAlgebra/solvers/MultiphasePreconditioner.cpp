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
 * @file PreconditionerMultiphase.hpp
 */

#include "MultiphasePreconditioner.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

#include <numeric>

namespace geosx
{

template< typename LAI >
MultiphasePreconditioner< LAI >::MultiphasePreconditioner( SchurApproximationType const schurType,
                                                           LinearSolverParameters parameters )
  : PreconditionerBase< LAI >(),
  m_schurType( schurType ),
  m_parameters( std::move( parameters ) )
{ }

template< typename LAI >
MultiphasePreconditioner< LAI >::~MultiphasePreconditioner()
{ }

template< typename LAI >
void MultiphasePreconditioner< LAI >::setPrimaryDofComponents( DofManager const & dofManager,
                                                               std::vector< DofManager::SubComponent > primaryDofs )
{
  m_dofComponentsPrimary = std::move( primaryDofs );
  m_dofComponentsOther = dofManager.filterDofs( m_dofComponentsPrimary );
}

namespace
{

/**
 * @brief Determine if two DOF subsets match, i.e. have the same support on the mesh.
 * @param dofManager
 * @param field1 name of first field
 * @param field2 name of second field
 * @return
 *
 * @note The detection procedure is not 100% reliable: it is possible in some rare cases for
 *       two dofs with disjoint supports but same location type to be identified as matching.
 *       An example would be a flow problem with N reservoir cells and N well segments.
 *
 * @todo Find a different way to detect this, maybe with user input
 */
bool checkMatchingDofs( DofManager const & dofManager,
                        string const & field1,
                        string const & field2 )
{
  GEOSX_ASSERT( dofManager.fieldExists( field1 ) );
  GEOSX_ASSERT( dofManager.fieldExists( field1 ) );

  bool result = field1 == field2;
#if 1
  result = result || ( dofManager.numGlobalSupport( field1 ) == dofManager.numGlobalSupport( field2 )
                       && dofManager.numLocalSupport( field1 ) == dofManager.numLocalSupport( field2 )
                       && dofManager.getLocation( field1 ) == dofManager.getLocation( field2 ) );
#endif
  return result;
}

}

template< typename LAI >
void MultiphasePreconditioner< LAI >::createReduction( DofManager const & dofManager )
{
  localIndex const numPrimaryFields = m_dofComponentsPrimary.size();
  localIndex const numOtherFields = m_dofComponentsOther.size();
  localIndex const numLocalPrimaryDof = m_prolongationOpPrimary.numLocalCols();
  localIndex const numLocalOtherDof = m_prolongationOpOther.numLocalCols();

  array1d< globalIndex > rowIndices;
  array1d< globalIndex > colIndices;
  array2d< real64 > values;

  localIndex const maxEntriesPerRow =
    std::accumulate( m_dofComponentsOther.begin(), m_dofComponentsOther.end(), 0,
                     []( localIndex const numComp, DofManager::SubComponent const & sc )
  { return std::max( numComp, sc.hiComp - sc.loComp ); } ) + 1;

  m_reductionOp.createWithLocalSize( numLocalPrimaryDof, numLocalPrimaryDof + numLocalOtherDof,
                                     maxEntriesPerRow, this->matrix().getComm() );
  m_reductionOp.open();

  globalIndex globalOffsetRow = m_reductionOp.ilower();
  for( localIndex blockRow = 0; blockRow < numPrimaryFields; ++blockRow )
  {
    DofManager::SubComponent const & dofPrim = m_dofComponentsPrimary[blockRow];
    localIndex const numLocalNodesPrim = dofManager.numLocalSupport( dofPrim.fieldName );
    localIndex const numCompPrim = dofPrim.hiComp - dofPrim.loComp;
    localIndex const numCompPrimField = dofManager.numComponents( dofPrim.fieldName );

    // 1. Assemble self-restriction block for the current primary dof

    rowIndices.resize( numCompPrim );
    colIndices.resize( numCompPrimField );
    values.resize( numCompPrim, numCompPrimField );
    values = 1.0;

    globalIndex globalOffsetCol = dofManager.globalOffset( dofPrim.fieldName );
    for( localIndex k = 0; k < numLocalNodesPrim; ++k )
    {
      for( localIndex c = 0; c < numCompPrim; ++c )
      {
        rowIndices[c] = globalOffsetRow + k * numCompPrim + c;
      }
      for( localIndex c = 0; c < numCompPrimField; ++c )
      {
        colIndices[c] = globalOffsetCol + k * numCompPrimField + c;
      }
      m_reductionOp.insert( rowIndices, colIndices, values );
    }

    // 2. Assemble reduction blocks for every matching "other" dof that is a different field

    for( localIndex blockCol = 0; blockCol < numOtherFields; ++blockCol )
    {
      DofManager::SubComponent const & dofOther = m_dofComponentsOther[blockCol];

      // Only dof pairs that have matching supports produce reduction operator entries
      if( checkMatchingDofs( dofManager, dofPrim.fieldName, dofOther.fieldName )
          && dofPrim.fieldName != dofOther.fieldName )
      {
        localIndex const numCompOther = dofOther.hiComp - dofOther.loComp;

        globalOffsetCol = dofManager.globalOffset( dofOther.fieldName );

        rowIndices.resize( numCompPrim );
        colIndices.resize( numCompOther );
        values.resize( numCompPrim, numCompOther );
        values = 1.0;

        for( localIndex k = 0; k < numLocalNodesPrim; ++k )
        {
          for( localIndex c = 0; c < numCompPrim; ++c )
          {
            rowIndices[c] = globalOffsetRow + k * numCompPrim + c;
          }
          for( localIndex c = 0; c < numCompOther; ++c )
          {
            colIndices[c] = globalOffsetCol + k * numCompOther + c;
          }

          m_reductionOp.insert( rowIndices, colIndices, values );
        }
      }
    }
    globalOffsetRow += numLocalNodesPrim * numCompPrim;
  }

  m_reductionOp.close();
}

namespace
{

// This constructs the block column summator: a tall matrix C  [ 1     ]
// s.t. A^T * C -> CS => a tall matrix where each nxn block    [   1   ]
// of rows is a diagonal block of A plus off-diagonal blocks   [     1 ]
// in the same column (both local and remote).                 [ 1     ]
//                                                             [   1   ]
// The matrix has the following shape:         ------->        [     1 ]
// (slightly more general for multiple dofs)                   [ ..... ]
template< typename VECTOR >
void createColsumOperator( DofManager const & dofManager,
                           std::vector< DofManager::SubComponent > const & dofComponents,
                           MPI_Comm const & comm,
                           array1d< VECTOR > & colsum )
{
  localIndex const numCompTotal =
    std::accumulate( dofComponents.begin(), dofComponents.end(), 0,
                     []( localIndex const numComp, DofManager::SubComponent const & sc )
  { return numComp + sc.hiComp - sc.loComp; } );

  localIndex const numLocalDof =
    std::accumulate( dofComponents.begin(), dofComponents.end(), 0,
                     [&]( localIndex const numDof, DofManager::SubComponent const & sc )
  { return numDof + dofManager.numLocalSupport( sc.fieldName ) * (sc.hiComp - sc.loComp); } );

  colsum.resize( numCompTotal );
  localIndex col = 0;

  localIndex rowOffset = 0;
  for( DofManager::SubComponent const & dof : dofComponents )
  {
    localIndex const numComp = dof.hiComp - dof.loComp;
    for( localIndex c = 0; c < numComp; ++c )
    {
      colsum[col + c].createWithLocalSize( numLocalDof, comm );
      colsum[col + c].zero();
      colsum[col + c].open();
    }

    localIndex const numLocalNodes = dofManager.numLocalSupport( dof.fieldName );
    globalIndex row = colsum[0].ilower();

    for( localIndex k = 0; k < numLocalNodes; ++k )
    {
      for( localIndex c = 0; c < numComp; ++c, ++row )
      {
        colsum[col + c].set( rowOffset + row, 1.0 );
      }
    }
    for( localIndex c = 0; c < numComp; ++c )
    {
      colsum[col + c].close();
    }
    col += numComp;
    rowOffset += numLocalNodes * numComp;
  }
}

}

template< typename LAI >
void MultiphasePreconditioner< LAI >::initialize( DofManager const & dofManager )
{
  MPI_Comm const & comm = this->matrix().getComm();

  dofManager.makeRestrictor( m_dofComponentsPrimary, comm, true, m_prolongationOpPrimary );
  dofManager.makeRestrictor( m_dofComponentsPrimary, comm, false, m_restrictionOpPrimary );
  dofManager.makeRestrictor( m_dofComponentsOther, comm, true, m_prolongationOpOther );
  dofManager.makeRestrictor( m_dofComponentsOther, comm, false, m_restrictionOpOther );

  localIndex const numLocalPrimaryDof = m_prolongationOpPrimary.numLocalCols();
  localIndex const numLocalOtherDof = m_prolongationOpOther.numLocalCols();
  GEOSX_LAI_ASSERT_EQ( numLocalPrimaryDof + numLocalOtherDof, this->matrix().numLocalCols() );

  localIndex const numCompPrimary =
    std::accumulate( m_dofComponentsPrimary.begin(), m_dofComponentsPrimary.end(), 0,
                     []( localIndex const numComp, DofManager::SubComponent const & sc )
  { return numComp + sc.hiComp - sc.loComp; } );

  localIndex const numCompOther =
    std::accumulate( m_dofComponentsOther.begin(), m_dofComponentsOther.end(), 0,
                     []( localIndex const numComp, DofManager::SubComponent const & sc )
  { return numComp + sc.hiComp - sc.loComp; } );

  // Primary system vectors
  m_rhsPrimary.createWithLocalSize( numLocalPrimaryDof, comm );
  m_solPrimary.createWithLocalSize( numLocalPrimaryDof, comm );

  // Generate sparsity pattern of Schur matrix
  createReduction( dofManager );

  // Colsum "operators" (multivectors)
  if( m_schurType == SchurApproximationType::COLSUM_BLOCK_DIAGONAL )
  {
    createColsumOperator( dofManager, m_dofComponentsPrimary, comm, m_colsumOpPrimary );
    createColsumOperator( dofManager, m_dofComponentsOther, comm, m_colsumOpOther );

    m_colsumResultPrimary.resize( numCompPrimary );
    for( localIndex c = 0; c < numCompPrimary; ++c )
    {
      m_colsumResultPrimary[c].createWithLocalSize( numLocalOtherDof, comm );
    }

    m_colsumResultOther.resize( numCompOther );
    for( localIndex c = 0; c < numCompOther; ++c )
    {
      m_colsumResultOther[c].createWithLocalSize( numLocalOtherDof, comm );
    }
  }

  // First stage preconditioner
  m_solverPrimary = LAI::createPreconditioner( m_parameters );
}

template< typename LAI >
void MultiphasePreconditioner< LAI >::computeColsumReduction( DofManager const & dofManager,
                                                              Matrix const & PS,
                                                              Matrix const & SS )
{
  localIndex const numPrimaryFields = m_dofComponentsPrimary.size();
  localIndex const numOtherFields   = m_dofComponentsOther.size();

  // Compute colsums and extract pointers to local data

  array1d< real64 const * > colsumPrimary( m_colsumOpPrimary.size() );
  for( localIndex j = 0; j < m_colsumOpPrimary.size(); ++j )
  {
    PS.applyTranspose( m_colsumOpPrimary[j], m_colsumResultPrimary[j] );
    colsumPrimary[j] = m_colsumResultPrimary[j].extractLocalVector();
  }

  array1d< real64 const * > colsumOther( m_colsumOpOther.size() );
  for( localIndex j = 0; j < m_colsumOpOther.size(); ++j )
  {
    SS.applyTranspose( m_colsumOpOther[j], m_colsumResultOther[j] );
    colsumOther[j] = m_colsumResultOther[j].extractLocalVector();
  }

  // 2. Compute PS * SS^{-1} block-by-block and scatter into reduction operator

  array2d< real64 > ps;
  array2d< real64 > ss;
  array2d< real64 > ss_inv;

  array1d< globalIndex > rowIndices;
  array1d< globalIndex > colIndices;
  array2d< real64 > values;

  localIndex compOffsetPrimary = 0;
  globalIndex globalRowOffset = m_reductionOp.ilower();

  m_reductionOp.zero();
  m_reductionOp.open();

  for( localIndex blockRow = 0; blockRow < numPrimaryFields; ++blockRow )
  {
    DofManager::SubComponent const & dofPrim = m_dofComponentsPrimary[blockRow];
    globalIndex globalColOffset = dofManager.globalOffset( dofPrim.fieldName );
    localIndex const numCompPrimary = dofPrim.hiComp - dofPrim.loComp;
    localIndex const numLocalSupport = dofManager.numLocalSupport( dofPrim.fieldName );
    localIndex const numCompPrimaryField = dofManager.numComponents( dofPrim.fieldName );

    for( localIndex k = 0; k < numLocalSupport; ++k )
    {
      globalIndex const rowIndex = globalRowOffset + k * numCompPrimary;
      globalIndex const colIndex = globalColOffset + k * numCompPrimaryField;
      for( localIndex c = 0; c < numCompPrimary; ++c )
      {
        m_reductionOp.set( rowIndex + c, colIndex + c, 1.0 );
      }
    }

    localIndex compOffsetOther = 0;
    localIndex rowOffsetOther = 0;
    for( localIndex blockCol = 0; blockCol < numOtherFields; ++blockCol )
    {
      DofManager::SubComponent const & dofOther = m_dofComponentsOther[blockCol];
      localIndex const numCompOther = dofOther.hiComp - dofOther.loComp;

      // Only dof pairs that have matching supports produce reduction operator entries
      if( checkMatchingDofs( dofManager, dofPrim.fieldName, dofOther.fieldName ) )
      {
        ps.resize( numCompPrimary, numCompOther );
        ss.resize( numCompOther, numCompOther );
        ss_inv.resize( numCompOther, numCompOther );
        values.resize( numCompPrimary, numCompOther );
        rowIndices.resize( numCompPrimary );
        colIndices.resize( numCompOther );

        // column indices are computed in a way that accounts for both
        // same-field and split-field formulations
        globalColOffset = dofManager.globalOffset( dofOther.fieldName );
        localIndex const colStride = numCompOther + ( dofPrim.fieldName == dofOther.fieldName ? numCompPrimary : 0 );
        localIndex const colShift = ( dofPrim.fieldName == dofOther.fieldName ? numCompPrimary : 0 );

        for( localIndex k = 0; k < numLocalSupport; ++k )
        {
          // extract and transpose off-diagonal block from PS
          for( localIndex j = 0; j < numCompPrimary; ++j )
          {
            for( localIndex i = 0; i < numCompOther; ++i )
            {
              ps( j, i ) = colsumPrimary[compOffsetPrimary + j][rowOffsetOther + k * numCompOther + i];
            }
          }

          // extract and transpose diagonal block from SS
          for( localIndex j = 0; j < numCompOther; ++j )
          {
            for( localIndex i = 0; i < numCompOther; ++i )
            {
              ss( j, i ) = colsumOther[compOffsetOther + j][rowOffsetOther + k * numCompOther + i];
            }
          }

          // Compute the inverse and product
          BlasLapackLA::matrixInverse( ss, ss_inv );
          BlasLapackLA::matrixMatrixMultiply( ps, ss_inv, values, -1.0 );

          // Scatter into global reduction operator
          for( localIndex i = 0; i < numCompPrimary; ++i )
          {
            rowIndices[i] = globalRowOffset + k * numCompPrimary + i;
          }
          for( localIndex j = 0; j < numCompOther; ++j )
          {
            colIndices[j] = globalColOffset + k * colStride + colShift + j;
          }

          m_reductionOp.set( rowIndices, colIndices, values );
        }
      }
      compOffsetOther += numCompOther;
      rowOffsetOther += dofManager.numLocalSupport( dofOther.fieldName ) * numCompOther;
    }
    compOffsetPrimary += numCompPrimary;
    globalRowOffset += numLocalSupport * numCompPrimary;
  }
  m_reductionOp.close();
}

template< typename LAI >
void MultiphasePreconditioner< LAI >::compute( Matrix const & mat,
                                               DofManager const & dofManager )
{
  bool const newSize = !this->ready() ||
                       mat.numGlobalRows() != this->numGlobalRows() ||
                       mat.numGlobalCols() != this->numGlobalRows();

  Base::compute( mat );

  // 1. Create reusable data structures (matrices/vectors)
  if( newSize )
  {
    initialize( dofManager );
  }

  // 2. Extract PS/SS blocks needed for reduction
  Matrix PS;
  mat.multiplyRAP( m_restrictionOpPrimary, m_prolongationOpOther, PS );
  Matrix SS;
  mat.multiplyPtAP( m_prolongationOpOther, SS );

  // 3. Construct reduction operator
  if( m_schurType == SchurApproximationType::COLSUM_BLOCK_DIAGONAL )
  {
    computeColsumReduction( dofManager, PS, SS );
  }
  else
  {
    GEOSX_ERROR( "Block-diagonal reduction not implemented" );
  }

  // 4. Construct reduced pressure matrix
  mat.multiplyRAP( m_reductionOp, m_prolongationOpPrimary, m_systemMatrixPrimary );

  // 5. Output matrices
  if( m_parameters.logLevel >= 3 )
  {
    //mat.write( "A.txt", LAIOutputFormat::MATRIX_MARKET );
    m_systemMatrixPrimary.write( "Ap.txt", LAIOutputFormat::MATRIX_MARKET );
    m_reductionOp.write( "Rp.txt", LAIOutputFormat::MATRIX_MARKET );
  }

  // 6. Setup pressure solver
  m_solverPrimary->compute( m_systemMatrixPrimary );
}

template< typename LAI >
void MultiphasePreconditioner< LAI >::apply( Vector const & src,
                                             Vector & dst ) const
{
  GEOSX_LAI_ASSERT( this->ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  // Apply schur complement restriction to obtain pressure system
  m_reductionOp.apply( src, m_rhsPrimary );

  // Apply the pressure solver
  m_solverPrimary->apply( m_rhsPrimary, m_solPrimary );

  // Prolongate pressure solution to full, using dst as accumulator
  m_prolongationOpPrimary.apply( m_solPrimary, dst );
}

template< typename LAI >
void MultiphasePreconditioner< LAI >::clear()
{
  Base::clear();
  m_solverPrimary->clear();
  // TODO: destroy all LA objects too
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MultiphasePreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MultiphasePreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MultiphasePreconditioner< PetscInterface >;
#endif

} //namespace geosx
