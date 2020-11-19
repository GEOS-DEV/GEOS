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
 * @file EpetraSuperLU_Dist.cpp
 */

#include "EpetraSuperLU_Dist.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/direct/Arnoldi.hpp"

#include <Epetra_FECrsMatrix.h>

namespace geosx
{

void EpetraConvertToSuperMatrix( EpetraMatrix const & matrix,
                                 SuperLU_Dist & SLUDData )
{
  int * offset;
  int * indices;
  real64 * values;
  GEOSX_LAI_CHECK_ERROR( matrix.unwrapped().ExtractCrsDataPointers( offset, indices, values ) );

  // contains the global ID of local columns
  globalIndex const * globalColumns = matrix.unwrapped().RowMatrixColMap().MyGlobalElements64();

  globalIndex const numGlobalRows = matrix.numGlobalRows();
  localIndex const numLocalRows = matrix.numLocalRows();
  localIndex const numLocalNonzeros = matrix.numLocalNonzeros();

  SLUDData.setNumGlobalRows( LvArray::integerConversion< int_t >( numGlobalRows ) );
  SLUDData.setNumGlobalCols( LvArray::integerConversion< int_t >( matrix.numGlobalCols() ) );

  SLUDData.resize( numLocalRows, numLocalNonzeros );
  for( localIndex i = 0; i <= numLocalRows; ++i )
  {
    SLUDData.rowPtr()[i] = LvArray::integerConversion< int_t >( offset[i] );
  }
  for( localIndex i = 0; i < numLocalNonzeros; ++i )
  {
    SLUDData.colIndices()[i] = LvArray::integerConversion< int_t >( globalColumns[indices[i]] );
    SLUDData.values()[i] = values[i];
  }

  SLUDData.createSuperMatrix( matrix.ilower() );
  SLUDData.setComm( matrix.getComm() );
}

namespace
{
class InverseNormalOperator : public LinearOperator< EpetraVector >
{
public:

  void set( EpetraMatrix const & matrix, SuperLU_Dist & SLUDData )
  {
    m_SLUDData = &SLUDData;
    m_comm = SLUDData.getComm();

    matrix.transpose( m_transposeMatrix );
    m_SLUDDataTransp.create( SLUDData.getParameters() );
    EpetraConvertToSuperMatrix( m_transposeMatrix, m_SLUDDataTransp );
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

  void apply( EpetraVector const & x, EpetraVector & y ) const override
  {
    m_SLUDData->solve( x.extractLocalVector(), y.extractLocalVector() );
    m_SLUDDataTransp.solve( y.extractLocalVector(), y.extractLocalVector() );
  }

private:

  SuperLU_Dist * m_SLUDData;

  MPI_Comm m_comm;

  EpetraMatrix m_transposeMatrix;

  mutable SuperLU_Dist m_SLUDDataTransp;
};
}

real64 EpetraSuperLU_DistCond( EpetraMatrix const & matrix, SuperLU_Dist & SLUDData )
{
  localIndex const numIterations = 4;

  using NormalOperator = NormalOperator< EpetraMatrix, EpetraVector >;
  NormalOperator normalOperator;
  normalOperator.set( matrix, matrix.getComm() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( normalOperator, numIterations );

  InverseNormalOperator inverseNormalOperator;
  inverseNormalOperator.set( matrix, SLUDData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseNormalOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
