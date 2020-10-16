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

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_EFEMKERNELSHELPER_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_EFEMKERNELSHELPER_HPP_

namespace geosx
{

namespace EFEMKernelsHelper
{

template<int NUM_NODES>
void ComputeHeavisideFunction( integer (& heaviside)[NUM_NODES],
		                       real64 (& X)[NUM_NODES][3],
							   arraySlice1d<real64 const> const & normalVector,
							   arraySlice1d<real64 const> const & elementCenter)
{
	for (int a=0; a < NUM_NODES; a++)
	{
		real64 distanceVector[3];
		LvArray::tensorOps::copy< 3 >( distanceVector, X[a] );
		LvArray::tensorOps::subtract< 3 >( distanceVector, elementCenter );

		heaviside[a] = LvArray::tensorOps::AiBi< 3 >( distanceVector, normalVector ) > 0 ? 1 : 0;
	}

}



template<int I_SIZE,
         int J_SIZE,
		 int NUM_NODES>
void AssembleStrainOperator( real64 ( & strainMatrix )[I_SIZE][J_SIZE],
		                     real64 ( & dNdX )[NUM_NODES][3]  )
{
  GEOSX_MARK_FUNCTION;
  LvArray::tensorOps::fill<I_SIZE, J_SIZE>(strainMatrix, 0);  //make 0
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
void AssembleCompatibilityOperator( real64 ( & compMatrix )[6][3],
		                            real64 ( & nVec )[3],
									real64 ( & tVec1 )[3],
									real64 ( & tVec2 )[3],
									real64 ( & heavisideFun )[NUM_NODES],
									real64 ( & dNdX )[NUM_NODES][3] )
{
  GEOSX_MARK_FUNCTION;

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point
  real64 mVec[3], heavisideFun;
  LvArray::tensorOps::fill<3>(mVec, 0);
  for( integer a=0; a<NUM_NODES; ++a )
  {
	LvArray::tensorOps::subtract<3>(mVec, dNdX[a]);
	LvArray::tensorOps::scale<3>(mVec, heavisideFun[a]);
  }

  // 2. fill in the operator itself

  LvArray::tensorOps::fill<6, 3>( compMatrix, 0 );

  real64 nDmSym[3][3], t1DmSym[3][3], t2DmSym[3][3];

  // sym(n dyadic m)
  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( nDmSym, mVec, nVec );
  LvArray::tensorOps::Rij_add_AiBj< 3, 3 >( nDmSym, nVec, mVec );
  LvArray::tensorOps::scale< 3, 3 >( nDmSym, 0.5 );

  // sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( t1DmSym, mVec, tVec1 );
  LvArray::tensorOps::Rij_add_AiBj< 3, 3 >( t1DmSym, tVec1, mVec );
  LvArray::tensorOps::scale< 3, 3 >( t1DmSym, 0.5 );

  LvArray::tensorOps::Rij_eq_AiBj< 3, 3 >( t2DmSym, mVec, tVec2 );
  LvArray::tensorOps::Rij_add_AiBj< 3, 3 >( t2DmSym, tVec2, mVec );
  LvArray::tensorOps::scale< 3, 3 >( t2DmSym, 0.5 );

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
      compMatrix( VoigtIndex, 0 ) += nDmSym [i][j];
      compMatrix( VoigtIndex, 1 ) += t1DmSym[i][j];
      compMatrix( VoigtIndex, 2 ) += t2DmSym[i][j];
    }
  }
}

void AssembleEquilibriumOperator( real64 ( & eqMatrix )[3][6],
		                          arraySlice1d<real64 const> const &  nVec,
							      arraySlice1d<real64 const> const &  tVec1,
								  arraySlice1d<real64 const> const &  tVec2,
                                  real64 const hInv );

}
} // geosx


#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_EFEMKERNELSHELPER_HPP_ */
