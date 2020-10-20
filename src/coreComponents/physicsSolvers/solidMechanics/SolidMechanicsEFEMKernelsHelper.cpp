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
 * @file EFEMKernelsHelper.cpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELSHELPER_CPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELSHELPER_CPP_

#include "SolidMechanicsEFEMKernelsHelper.hpp"

namespace geosx
{

namespace SolidMechanicsEFEMKernelsHelper
{

void assembleEquilibriumOperator( real64 ( & eqMatrix )[3][6],
                                  arraySlice1d< real64 const > const & nVec,
                                  arraySlice1d< real64 const > const & tVec1,
                                  arraySlice1d< real64 const > const & tVec2,
                                  real64 const hInv )
{
  GEOSX_MARK_FUNCTION;

  LvArray::tensorOps::fill< 3, 6 >( eqMatrix, 0 );

  // Define the dyadic producs
  // nDn := nVec_dyadic_nVec,
  // t1DnSym:= sym(tVec1_dyadic_nVec)
  // t2DnSym:= sym(tVec2_dyadic_nVec)
  real64 nDn[3][3], t1DnSym[3][3], t2DnSym[3][3];

  // n dyadic n
  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( nDn, nVec, nVec );

  // sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( t1DnSym, nVec, tVec1 );
  LvArray::tensorOps::Rij_add_AiBj< 3, 3 >( t1DnSym, tVec1, nVec );
  LvArray::tensorOps::scale< 3, 3 >( t1DnSym, 0.5 );

  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( t2DnSym, nVec, tVec2 );
  LvArray::tensorOps::Rij_add_AiBj< 3, 3 >( t2DnSym, tVec2, nVec );
  LvArray::tensorOps::scale< 3, 3 >( t2DnSym, 0.5 );

  int VoigtIndex;

  for( int i=0; i < 3; ++i )
  {
    for( int j=0; j < 3; ++j )
    {
      if( i == j )
      {
        VoigtIndex = 1;
      }
      else
      {
        VoigtIndex = 6 - i - j;
      }
      eqMatrix[0][VoigtIndex] += nDn     [i][j];
      eqMatrix[1][VoigtIndex] += t1DnSym [i][j];
      eqMatrix[2][VoigtIndex] += t2DnSym [i][j];
    }
  }
  LvArray::tensorOps::scale< 3, 6 >( eqMatrix, -hInv );
}


void computeTraction( real64 ( & dispJump )[3],
		              real64 contactCoeff,
					  real64 ( & tractionVector )[3],
					  real64 ( & dTractiondw )[3][3] )
{
// check if fracture is open
  bool open = dispJump[0] >= 0 ? true : false;

  LvArray::tensorOps::fill<3>(tractionVector, 0);
  LvArray::tensorOps::fill<3, 3>(dTractiondw, 0);

  if( open )
  {
    tractionVector[0] = 1e5;
  }
  else
  {
    // Contact through penalty condition.
    tractionVector[0] = contactCoeff * dispJump[0];
    dTractiondw[0][0] = contactCoeff;
  }

}

} // namespace EFEMKernelsHelper
} // namespace geosx



#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_EFEMKERNELSHELPER_CPP_ */
