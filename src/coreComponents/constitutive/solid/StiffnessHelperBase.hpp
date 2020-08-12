/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_STIFFNESSHELPERBASE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_STIFFNESSHELPERBASE_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{
namespace constitutive
{


struct StiffnessHelperBase
{
  template< int NUM_SUPPORT_POINTS,
            int NUM_DOF_PER_SUPPORT_POINT,
            typename BASIS_GRADIENT,
            typename CBF >
  void BTDB( BASIS_GRADIENT const & dNdX,
             real64 ( &localJacobian )[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT],
             CBF && callbackFunction );

  template< int NUM_SUPPORT_POINTS,
            int NUM_DOF_PER_SUPPORT_POINT,
            typename BASIS_GRADIENT,
            typename CBF >
  void diagBTDB( BASIS_GRADIENT const & dNdX,
                 real64 ( &diagLocalJacobian )[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT],
                 CBF && callbackFunction );

};


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT,
          typename CBF >
void StiffnessHelperBase::BTDB( BASIS_GRADIENT const & dNdX,
                                real64 (& localJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT],
                                CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const dNdXadNdXb[3][3] = { { dNdX[a][0] * dNdX[b][0], dNdX[a][0] * dNdX[b][1], dNdX[a][0] * dNdX[b][2] },
        { dNdX[a][1] * dNdX[b][0], dNdX[a][1] * dNdX[b][1], dNdX[a][1] * dNdX[b][2] },
        { dNdX[a][2] * dNdX[b][0], dNdX[a][2] * dNdX[b][1], dNdX[a][2] * dNdX[b][2] } };
      callbackFunction( a, b, dNdXadNdXb, localJacobian );
    }
  }

  for( int row=0; row<NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT; ++row )
  {
    for( int col=0; col<row; ++col )
    {
      localJacobian[row][col] = localJacobian[col][row];
    }
  }
}


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT,
          typename CBF >
void StiffnessHelperBase::diagBTDB( BASIS_GRADIENT const & dNdX,
                                    real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT],
                                    CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const dNdXdNdX[3] = { dNdX[a][0] * dNdX[a][0], dNdX[a][1] * dNdX[a][1], dNdX[a][2] * dNdX[a][2] };

    callbackFunction( a, dNdXdNdX, diagLocalJacobian );
  }
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_STIFFNESSHELPERBASE_HPP_ */
