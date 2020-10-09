/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreSuperLU_Dist.cpp
 */

#include "HypreSuperLU_Dist.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/direct/Arnoldi.hpp"

#include "_hypre_parcsr_mv.h"

namespace geosx
{

void HypreConvertToSuperMatrix( HypreMatrix const & matrix,
                                hypre_CSRMatrix * & localMatrix,
                                SuperLU_Dist & SLUDData )
{
  // Merge diag and offd into one matrix (global ids)
  localMatrix = hypre_MergeDiagAndOffd( matrix.unwrapped() );

  globalIndex const numGlobalRows = matrix.numGlobalRows();
  localIndex const numLocalRows = matrix.numLocalRows();
  localIndex const numLocalNonzeros = matrix.numLocalNonzeros();

  SLUDData.setNumGlobalRows( LvArray::integerConversion< int_t >( numGlobalRows ) );
  SLUDData.setNumLocalRows( LvArray::integerConversion< int_t >( numLocalRows ) );

  HYPRE_Int const * const hypreI = hypre_CSRMatrixI( localMatrix );
  SLUDData.createRowPtr( numLocalRows );
  for( localIndex i = 0; i <= numLocalRows; ++i )
  {
    SLUDData.rowPtr()[i] = LvArray::integerConversion< int_t >( hypreI[i] );
  }
  SLUDData.setColIndices( toSuperLU_intT( hypre_CSRMatrixBigJ( localMatrix ) ) );
  SLUDData.setValues( hypre_CSRMatrixData( localMatrix ) );

  dCreate_CompRowLoc_Matrix_dist( &SLUDData.mat(),
                                  toSuperLU_intT( numGlobalRows ),
                                  toSuperLU_intT( numGlobalRows ),
                                  toSuperLU_intT( numLocalNonzeros ),
                                  toSuperLU_intT( numLocalRows ),
                                  toSuperLU_intT( matrix.ilower() ),
                                  SLUDData.values(),
                                  SLUDData.colIndices(),
                                  SLUDData.rowPtr(),
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );

  SLUDData.setComm( matrix.getComm() );
}

void HypreDestroyAdditionalData( hypre_CSRMatrix * & localMatrix )
{
  // From HYPRE SuperLU_Dist interfaces (superlu.c)
  // SuperLU frees assigned data, so set them to null before
  // calling hypre_CSRMatrixdestroy on localMatrix to avoid memory errors.
  hypre_CSRMatrixData( localMatrix ) = NULL;
  hypre_CSRMatrixBigJ( localMatrix ) = NULL;
  hypre_CSRMatrixDestroy( localMatrix );
}

namespace
{
class InverseOperator
{
public:

  ~InverseOperator()
  {
    HypreDestroyAdditionalData( m_localMatrix );
  }

  void set( HypreMatrix const & matrix, SuperLU_Dist & SLUDData )
  {
    m_SLUDData = &SLUDData;
    m_comm = SLUDData.getComm();

    matrix.transpose( m_transposeMatrix );
    m_SLUDDataTransp.create( SLUDData.getParameters() );
    HypreConvertToSuperMatrix( m_transposeMatrix, m_localMatrix, m_SLUDDataTransp );
    m_SLUDDataTransp.setup();
  }

  globalIndex globalSize() const
  {
    return LvArray::integerConversion< globalIndex >( m_SLUDData->numGlobalRows() );
  }

  localIndex localSize() const
  {
    return LvArray::integerConversion< localIndex >( m_SLUDData->numLocalRows() );
  }

  MPI_Comm const & getComm() const
  {
    return m_comm;
  }

  void apply( HypreVector const & x, HypreVector & y ) const
  {
    m_SLUDData->solve( x.extractLocalVector(), y.extractLocalVector() );
    m_SLUDDataTransp.solve( y.extractLocalVector(), y.extractLocalVector() );
  }

private:

  SuperLU_Dist * m_SLUDData;

  MPI_Comm m_comm;

  HypreMatrix m_transposeMatrix;

  hypre_CSRMatrix * m_localMatrix;

  mutable SuperLU_Dist m_SLUDDataTransp;
};
}

real64 HypreSuperLU_DistCond( HypreMatrix const & matrix, SuperLU_Dist & SLUDData )
{
  localIndex const numIterations = 4;

  using DirectOperator = DirectOperator< HypreMatrix, HypreVector >;
  DirectOperator directOperator;
  directOperator.set( matrix );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue< HypreVector, DirectOperator >( directOperator, numIterations );

  InverseOperator inverseOperator;
  inverseOperator.set( matrix, SLUDData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue< HypreVector, InverseOperator >( inverseOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
