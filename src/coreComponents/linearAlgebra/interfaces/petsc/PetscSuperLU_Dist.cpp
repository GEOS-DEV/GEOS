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
 * @file PetscSuperLU_Dist.cpp
 */

#include "PetscSuperLU_Dist.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/direct/Arnoldi.hpp"

namespace geosx
{

void PetscConvertToSuperMatrix( PetscMatrix const & matrix,
                                Mat & localMatrix,
                                SuperLU_Dist & SLUDData )
{
  GEOSX_LAI_CHECK_ERROR( MatMPIAIJGetLocalMat( matrix.unwrapped(), MAT_INITIAL_MATRIX, &localMatrix ) );

  PetscInt numRows;
  const PetscInt * ia;
  const PetscInt * ja;
  PetscBool status;
  GEOSX_LAI_CHECK_ERROR( MatGetRowIJ( localMatrix, 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
  if( !status )
  {
    GEOSX_ERROR( "ConvertPetscToSuperMatrix: MatGetRowIJ reported an error." );
  }

  real64 * array;
  GEOSX_LAI_CHECK_ERROR( MatSeqAIJGetArray( localMatrix, &array ) );

  globalIndex const numGlobalRows = matrix.numGlobalRows();
  localIndex const numLocalRows = matrix.numLocalRows();
  localIndex const numLocalNonzeros = matrix.numLocalNonzeros();

  SLUDData.setNumGlobalRows( LvArray::integerConversion< int_t >( numGlobalRows ) );
  SLUDData.setNumLocalRows( LvArray::integerConversion< int_t >( numLocalRows ) );

  SLUDData.createRowPtr( numLocalRows );
  for( localIndex i = 0; i <= numLocalRows; ++i )
  {
    SLUDData.rowPtr()[i] = LvArray::integerConversion< int_t >( ia[i] );
  }
  SLUDData.createColIndices( numLocalNonzeros );
  SLUDData.createValues( numLocalNonzeros );
  for( localIndex i = 0; i < numLocalNonzeros; ++i )
  {
    SLUDData.colIndices()[i] = LvArray::integerConversion< int_t >( ja[i] );
    SLUDData.values()[i] = array[i];
  }
  GEOSX_LAI_CHECK_ERROR( MatRestoreRowIJ( localMatrix, 0, PETSC_FALSE, PETSC_FALSE, &numRows, &ia, &ja, &status ) );
  if( !status )
  {
    GEOSX_ERROR( "ConvertPetscToSuperMatrix: MatGetRowIJ reported an error." );
  }
  GEOSX_LAI_CHECK_ERROR( MatSeqAIJRestoreArray( localMatrix, &array ) );

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

void PetscDestroyAdditionalData( Mat & localMatrix )
{
  GEOSX_LAI_CHECK_ERROR( MatDestroy( &localMatrix ) );
}

namespace
{
class InverseOperator
{
public:

  ~InverseOperator()
  {
    PetscDestroyAdditionalData( m_localMatrix );
  }

  void set( PetscMatrix const & matrix, SuperLU_Dist & SLUDData )
  {
    m_SLUDData = &SLUDData;
    m_comm = SLUDData.getComm();

    matrix.transpose( m_transposeMatrix );
    m_SLUDDataTransp.create( SLUDData.getParameters() );
    PetscConvertToSuperMatrix( m_transposeMatrix, m_localMatrix, m_SLUDDataTransp );
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

  void apply( PetscVector const & x, PetscVector & y ) const
  {
    m_SLUDData->solve( x.extractLocalVector(), y.extractLocalVector() );
    m_SLUDDataTransp.solve( y.extractLocalVector(), y.extractLocalVector() );
  }

private:

  SuperLU_Dist * m_SLUDData;

  MPI_Comm m_comm;

  PetscMatrix m_transposeMatrix;

  Mat m_localMatrix;

  mutable SuperLU_Dist m_SLUDDataTransp;
};
}

real64 PetscSuperLU_DistCond( PetscMatrix const & matrix, SuperLU_Dist & SLUDData )
{
  localIndex const numIterations = 4;

  using DirectOperator = DirectOperator< PetscMatrix, PetscVector >;
  DirectOperator directOperator;
  directOperator.set( matrix );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue< PetscVector, DirectOperator >( directOperator, numIterations );

  InverseOperator inverseOperator;
  inverseOperator.set( matrix, SLUDData );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue< PetscVector, InverseOperator >( inverseOperator, numIterations );

  return sqrt( lambdaDirect * lambdaInverse );
}

}
