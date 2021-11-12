/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscExport.cpp
 */

#include "PetscExport.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"

#include <petsc.h>

namespace geosx
{

PetscExport::PetscExport() = default;

PetscExport::PetscExport( PetscMatrix const & mat,
                          integer targetRank )
  : m_targetRank( targetRank )
{
  globalIndex const numGlobalRows = mat.numGlobalRows();
  localIndex const numLocalRows = mat.numLocalRows();

  int const rank = MpiWrapper::commRank( mat.getComm() );
  localIndex const N = ( rank == m_targetRank ) ? numGlobalRows : 0;

  // create vector scatter context
  Vec tmpGlobal;
  Vec tmpLocal;
  GEOSX_LAI_CHECK_ERROR( VecCreateMPI( mat.getComm(), numLocalRows, numGlobalRows, &tmpGlobal ) );
  GEOSX_LAI_CHECK_ERROR( VecCreateSeq( PETSC_COMM_SELF, N, &tmpLocal ) );
  GEOSX_LAI_CHECK_ERROR( ISCreateStride( PETSC_COMM_SELF, N, 0, 1, &m_indexSet ) );
  GEOSX_LAI_CHECK_ERROR( VecScatterCreate( tmpGlobal, m_indexSet, tmpLocal, m_indexSet, &m_scatter ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &tmpGlobal ) );
  GEOSX_LAI_CHECK_ERROR( VecDestroy( &tmpLocal ) );
}

PetscExport::~PetscExport()
{
  ISDestroy( &m_indexSet );
  VecScatterDestroy( &m_scatter );
}

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void PetscExport::exportCRS( PetscMatrix const & mat,
                             arrayView1d< OFFSET_TYPE > const & rowOffsets,
                             arrayView1d< COLUMN_TYPE > const & colIndices,
                             arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( mat.getComm() );

  // import on target rank if needed, or extract diag+offdiag part in each rank
  Mat * submat; // needed by MatCreateSubMatrices API
  Mat localMatrix;

  rowOffsets.move( LvArray::MemorySpace::host, false );
  colIndices.move( LvArray::MemorySpace::host, false );
  values.move( LvArray::MemorySpace::host, false );

  if( m_targetRank < 0 )
  {
    GEOSX_LAI_CHECK_ERROR( MatMPIAIJGetLocalMat( mat.unwrapped(), MAT_INITIAL_MATRIX, &localMatrix ) );
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( MatCreateSubMatrices( mat.unwrapped(), rank == m_targetRank ? 1 : 0,
                                                 &m_indexSet, &m_indexSet, MAT_INITIAL_MATRIX, &submat ) );
    localMatrix = rank == m_targetRank ? submat[0] : nullptr;
  }

  // export the raw CRS data
  if( m_targetRank < 0 || m_targetRank == rank )
  {
    PetscInt numRows;
    PetscInt const * ia;
    PetscInt const * ja;
    PetscBool status;
    GEOSX_LAI_CHECK_ERROR( MatGetRowIJ( localMatrix, 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
    GEOSX_ERROR_IF( !status, "PetscExtract: MatGetRowIJ reported an error" );

    real64 const * va;
    GEOSX_LAI_CHECK_ERROR( MatSeqAIJGetArrayRead( localMatrix, &va ) );

    MatInfo info;
    GEOSX_LAI_CHECK_ERROR( MatGetInfo( localMatrix, MAT_LOCAL, &info ) );
    PetscInt const numNz = static_cast< PetscInt >( info.nz_used );

    std::transform( ia, ia + numRows + 1, rowOffsets.data(), LvArray::integerConversion< OFFSET_TYPE, PetscInt > );
    std::transform( ja, ja + numNz, colIndices.data(), LvArray::integerConversion< COLUMN_TYPE, PetscInt > );
    std::copy( va, va + numNz, values.data() );

    GEOSX_LAI_CHECK_ERROR( MatSeqAIJRestoreArrayRead( localMatrix, &va ) );
    GEOSX_LAI_CHECK_ERROR( MatRestoreRowIJ( localMatrix, 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
  }

  if( m_targetRank < 0 )
  {
    GEOSX_LAI_CHECK_ERROR( MatDestroy( &localMatrix ));
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( MatDestroySubMatrices( rank == m_targetRank ? 1 : 0, &submat ) );
  }
}

void PetscExport::exportVector( PetscVector const & vec,
                                arrayView1d< real64 > const & values ) const
{
  values.move( LvArray::MemorySpace::host, false );
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    localIndex const N = ( rank == m_targetRank ) ? vec.globalSize() : 0;

    Vec localVector;
    GEOSX_LAI_CHECK_ERROR( VecCreateSeqWithArray( PETSC_COMM_SELF, 1, N, values.data(), &localVector ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterBegin( m_scatter, vec.unwrapped(), localVector, INSERT_VALUES, SCATTER_FORWARD ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterEnd( m_scatter, vec.unwrapped(), localVector, INSERT_VALUES, SCATTER_FORWARD ) );
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &localVector ) );
  }
  else
  {
    arrayView1d< real64 const > const data = vec.values();
    data.move( LvArray::MemorySpace::host, false );
    std::copy( data.begin(), data.end(), values.data() );
  }
}

void PetscExport::importVector( arrayView1d< const real64 > const & values,
                                PetscVector & vec ) const
{
  values.move( LvArray::MemorySpace::host, false );
  if( m_targetRank >= 0 )
  {
    int const rank = MpiWrapper::commRank( vec.getComm() );
    localIndex const N = ( rank == m_targetRank ) ? vec.globalSize() : 0;

    // Note: PETSc creates a vector wrapper around values by taking a pointer-to-const and casting away const-ness (!)
    //       This would lead to UB if we did anything that modified values, but here we only scatter (read) them.
    Vec localVector;
    GEOSX_LAI_CHECK_ERROR( VecCreateSeqWithArray( PETSC_COMM_SELF, 1, N, values.data(), &localVector ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterBegin( m_scatter, localVector, vec.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterEnd( m_scatter, localVector, vec.unwrapped(), INSERT_VALUES, SCATTER_REVERSE ) );
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &localVector ) );
  }
  else
  {
    arrayView1d< real64 > const data = vec.open();
    std::copy( values.data(), values.data() + vec.localSize(), data.begin() );
    vec.close();
  }
}

/**
 * Explicit template instantiation macro for PetscExport::exportCRS.
 * We need to explicitly instantiate this function template because:
 * - we want CRS consumers to specify their own destination buffer types;
 * - we're "hiding" PETSc headers from being included by consumer code.
 */
#define INST_PETSC_EXPORT_CRS( OFFSET_TYPE, COLUMN_TYPE ) \
  template void \
  PetscExport::exportCRS< OFFSET_TYPE, COLUMN_TYPE >( PetscMatrix const &, \
                                                      arrayView1d< OFFSET_TYPE > const &, \
                                                      arrayView1d< COLUMN_TYPE > const &, \
                                                      arrayView1d< real64 > const & ) const

// Add other instantiations as needed (only use built-in types)
INST_PETSC_EXPORT_CRS( int, int );
INST_PETSC_EXPORT_CRS( long, long );
INST_PETSC_EXPORT_CRS( long long, long long );

#undef INST_PETSC_EXPORT_CRS

} // namespace geosx
