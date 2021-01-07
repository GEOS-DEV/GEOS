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
                                SuperLU_Dist & SLUDData )
{
  // Merge diag and offd into one matrix (global ids)
  hypre_CSRMatrix * localMatrix = hypre_MergeDiagAndOffd( matrix.unwrapped() );

  globalIndex const numGlobalRows = matrix.numGlobalRows();
  localIndex const numLocalRows = matrix.numLocalRows();
  localIndex const numLocalNonzeros = matrix.numLocalNonzeros();

  SLUDData.setNumGlobalRows( LvArray::integerConversion< int_t >( numGlobalRows ) );
  SLUDData.setNumGlobalCols( LvArray::integerConversion< int_t >( matrix.numGlobalCols() ) );

  HYPRE_Int const * const hypreI = hypre_CSRMatrixI( localMatrix );
  HYPRE_BigInt const * const hypreJ = hypre_CSRMatrixBigJ( localMatrix );
  real64 const * const hypreData = hypre_CSRMatrixData( localMatrix );

  SLUDData.resize( numLocalRows, numLocalNonzeros );
  for( localIndex i = 0; i <= numLocalRows; ++i )
  {
    SLUDData.rowPtr()[i] = LvArray::integerConversion< int_t >( hypreI[i] );
  }
  for( localIndex i = 0; i < numLocalNonzeros; ++i )
  {
    SLUDData.colIndices()[i] = LvArray::integerConversion< int_t >( hypreJ[i] );
    SLUDData.values()[i] = hypreData[i];
  }

  SLUDData.createSuperMatrix( matrix.ilower() );
  SLUDData.setComm( matrix.getComm() );

  // From HYPRE SuperLU_Dist interfaces (superlu.c)
  // SuperLU frees assigned data, so set them to null before
  // calling hypre_CSRMatrixdestroy on localMatrix to avoid memory errors.
  hypre_CSRMatrixData( localMatrix ) = NULL;
  hypre_CSRMatrixBigJ( localMatrix ) = NULL;
  hypre_CSRMatrixDestroy( localMatrix );
}

namespace
{
class InverseNormalOperator : public LinearOperator< HypreVector >
{
public:

  void set( HypreMatrix const & matrix, SuperLU_Dist & SLUDData )
  {
    m_SLUDData = &SLUDData;
    m_comm = SLUDData.getComm();

    matrix.transpose( m_transposeMatrix );
    m_SLUDDataTransp.create( SLUDData.getParameters() );
    HypreConvertToSuperMatrix( m_transposeMatrix, m_SLUDDataTransp );
    m_SLUDDataTransp.setup();
  }

  globalIndex numGlobalRows() const override
  {
    return LvArray::integerConversion< globalIndex >( m_SLUDData->numGlobalRows() );
  }

  globalIndex numGlobalCols() const override
  {
    return LvArray::integerConversion< globalIndex >( m_SLUDData->numGlobalCols() );
  }

  localIndex numLocalRows() const
  {
    return LvArray::integerConversion< localIndex >( m_SLUDData->numLocalRows() );
  }

  MPI_Comm const & getComm() const
  {
    return m_comm;
  }

  void apply( HypreVector const & x, HypreVector & y ) const override
  {
    m_SLUDData->solve( x.extractLocalVector(), y.extractLocalVector() );
    m_SLUDDataTransp.solve( y.extractLocalVector(), y.extractLocalVector() );
  }

private:

  SuperLU_Dist * m_SLUDData;

  MPI_Comm m_comm;

  HypreMatrix m_transposeMatrix;

  mutable SuperLU_Dist m_SLUDDataTransp;
};
}

real64 HypreSuperLU_DistCond( HypreMatrix const & matrix, SuperLU_Dist & SLUDData )
{
  localIndex const numIterations = 4;

  using NormalOperator = NormalOperator< HypreMatrix, HypreVector >;
  NormalOperator normalOperator;
  normalOperator.set( matrix, matrix.getComm() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( normalOperator, numIterations );

  InverseNormalOperator inverseNormalOperator;
  inverseNormalOperator.set( matrix, SLUDData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseNormalOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
