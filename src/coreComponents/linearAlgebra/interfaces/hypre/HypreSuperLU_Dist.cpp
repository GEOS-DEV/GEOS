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

}
