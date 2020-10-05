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

  SLUDData.createRowPtr( numLocalRows );
  for( localIndex i = 0; i <= numLocalRows; ++i )
  {
    SLUDData.rowPtr()[i] = LvArray::integerConversion< int_t >( offset[i] );
  }
  SLUDData.createColIndices( numLocalNonzeros );
  SLUDData.createValues( numLocalNonzeros );
  for( localIndex i = 0; i < numLocalNonzeros; ++i )
  {
    SLUDData.colIndices()[i] = LvArray::integerConversion< int_t >( globalColumns[indices[i]] );
    SLUDData.values()[i] = values[i];
  }

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

}
