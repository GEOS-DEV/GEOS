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
 * @file EFEMKernelsHelper.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELSHELPER_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELSHELPER_HPP_


#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/GeosxMacros.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace SolidMechanicsEFEMKernelsHelper
{

template< int NUM_NODES >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeHeavisideFunction( integer (& heaviside)[NUM_NODES],
                               real64 (& X)[NUM_NODES][3],
                               arraySlice1d< real64 const > const normalVector,
                               arraySlice1d< real64 const > const elementCenter )
{
  for( int a=0; a < NUM_NODES; a++ )
  {
    real64 distanceVector[3];
    LvArray::tensorOps::copy< 3 >( distanceVector, X[a] );
    LvArray::tensorOps::subtract< 3 >( distanceVector, elementCenter );

    heaviside[a] = LvArray::tensorOps::AiBi< 3 >( distanceVector, normalVector ) > 0 ? 1 : 0;
  }

}



template< int I_SIZE,
          int J_SIZE,
          int NUM_NODES >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void assembleStrainOperator( real64 ( & strainMatrix )[I_SIZE][J_SIZE],
                             real64 ( & dNdX )[NUM_NODES][3] )
{
  LvArray::tensorOps::fill< I_SIZE, J_SIZE >( strainMatrix, 0 );  //make 0
  for( int a=0; a < NUM_NODES; ++a )
  {

    strainMatrix[0][a*3 + 0] = dNdX[a][0];
    strainMatrix[1][a*3 + 1] = dNdX[a][1];
    strainMatrix[2][a*3 + 2] = dNdX[a][2];

    strainMatrix[3][a*3 + 1] = dNdX[a][2];
    strainMatrix[3][a*3 + 2] = dNdX[a][1];

    strainMatrix[4][a*3 + 0] = dNdX[a][2];
    strainMatrix[4][a*3 + 2] = dNdX[a][0];

    strainMatrix[5][a*3 + 0] = dNdX[a][1];
    strainMatrix[5][a*3 + 1] = dNdX[a][0];
  }
}

template< int NUM_NODES >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void assembleCompatibilityOperator( real64 ( & compMatrix )[6][3],
                                    arraySlice1d< real64 const > const & nVec,
                                    arraySlice1d< real64 const > const & tVec1,
                                    arraySlice1d< real64 const > const & tVec2,
                                    integer ( & heavisideFun )[NUM_NODES],
                                    real64 ( & dNdX )[NUM_NODES][3] )
{
  //GEOSX_MARK_FUNCTION;

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point

  real64 mVec[3];
  LvArray::tensorOps::fill< 3 >( mVec, 0 );
  for( integer a=0; a<NUM_NODES; ++a )
  {
    mVec[0] -= dNdX[a][0] * heavisideFun[a];
    mVec[1] -= dNdX[a][1] * heavisideFun[a];
    mVec[2] -= dNdX[a][2] * heavisideFun[a];
  }


  // 2. fill in the operator itself

  real64 nDmSym[6], t1DmSym[6], t2DmSym[6];

  // sym(n dyadic m), sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::Rij_eq_AiBj_plus_AjBi< 3 >( nDmSym, mVec, nVec );
  LvArray::tensorOps::Rij_eq_AiBj_plus_AjBi< 3 >( t1DmSym, mVec, tVec1 );
  LvArray::tensorOps::Rij_eq_AiBj_plus_AjBi< 3 >( t2DmSym, mVec, tVec2 );

  // if compMatrix was real64[3][6], then you could just pass
  // compMatrix[0], compMatrix[1], compMatrix[2] instead of creating
  // nDmSum, t1DmSum, t2DmSum.

  compMatrix[0][0] = 0.5 * nDmSym[0];
  compMatrix[1][0] = 0.5 * nDmSym[1];
  compMatrix[2][0] = 0.5 * nDmSym[2];
  compMatrix[3][0] = nDmSym[3];
  compMatrix[4][0] = nDmSym[4];
  compMatrix[5][0] = nDmSym[5];

  compMatrix[0][1] = 0.5 * t1DmSym[0];
  compMatrix[1][1] = 0.5 * t1DmSym[1];
  compMatrix[2][1] = 0.5 * t1DmSym[2];
  compMatrix[3][1] = t1DmSym[3];
  compMatrix[4][1] = t1DmSym[4];
  compMatrix[5][1] = t1DmSym[5];

  compMatrix[0][2] = 0.5 * t2DmSym[0];
  compMatrix[1][2] = 0.5 * t2DmSym[1];
  compMatrix[2][2] = 0.5 * t2DmSym[2];
  compMatrix[3][2] = t2DmSym[3];
  compMatrix[4][2] = t2DmSym[4];
  compMatrix[5][2] = t2DmSym[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void assembleEquilibriumOperator( real64 ( & eqMatrix )[3][6],
                                  arraySlice1d< real64 const > const nVec,
                                  arraySlice1d< real64 const > const tVec1,
                                  arraySlice1d< real64 const > const tVec2,
                                  real64 const hInv )
{
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

/*
 * @brief Computes traction and derivative on each fracture segment.
 * @param constitutiveManager constant pointer to the constitutive mamanger
 * @param dispJump displacement jump
 * @param tractionVector traction vector
 * @param dTdw Derivative of the traction w.r.t. the jump.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeTraction( real64 ( & dispJump )[3],
                      real64 contactCoeff,
                      real64 ( & tractionVector )[3],
                      real64 ( & dTractiondw )[3][3] )
{
  // check if fracture is open
  bool open = dispJump[0] >= 0 ? true : false;

  LvArray::tensorOps::fill< 3 >( tractionVector, 0 );
  LvArray::tensorOps::fill< 3, 3 >( dTractiondw, 0 );

  if( open )
  {
    tractionVector[0] = 1.0e5;
  }
  else
  {
    // Contact through penalty condition.
    tractionVector[0] = contactCoeff * dispJump[0];
    dTractiondw[0][0] = contactCoeff;
  }
}

}
} // geosx


#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_EFEMKERNELSHELPER_HPP_ */
