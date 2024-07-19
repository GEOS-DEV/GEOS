/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TaperKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERKERNEL_HPP_

#include "WaveSolverUtils.hpp"

namespace geos
{

using wsCoordType = WaveSolverUtils::wsCoordType;

struct TaperKernel
{

  /**
   * @brief Compute coefficients for the taper layers
   * @tparam EXEC_POLICY the execution policy
   * @param[in] size the number of nodes
   * @param[in] xMin coordinate limits of the inner taper boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner taper boundaries, right-back-bottom
   * @param[in] dMin Taper thickness, left-front-top
   * @param[in] dMax Taper thickness, right-back-bottom
   * @param[in] dt time-step
   * @param[in] vMax
   * @param[in] r desired reflectivity of the Taper
   * @param[out] taperCoeff array which contains the taper coefficient on each node (which will be equal to 1 when we are outside of the
   * taper layers)
   */
  template< typename EXEC_POLICY >
  static void
  computeTaperCoeff( localIndex const size,
                     arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     R1Tensor32 const xMin,
                     R1Tensor32 const xMax,
                     R1Tensor32 const dMin,
                     R1Tensor32 const dMax,
                     real32 const dt,
                     real32 const vMax,
                     real32 const r,
                     arrayView1d< real32 > const taperCoeff )
  {
    /// Loop over elements in the subregion, 'l' is the element index within the target set

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      real32 const dxMin = -((3*vMax)/(2*dMin[0]))*log( r )*pow((dMin[0]-nodeCoords[a][0])/dMin[0], 2 );
      real32 const dyMin = -((3*vMax)/(2*dMin[1]))*log( r )*pow((dMin[1]-nodeCoords[a][1])/dMin[1], 2 );
      real32 const dzMin = -((3*vMax)/(2*dMin[2]))*log( r )*pow((dMin[2]-nodeCoords[a][2])/dMin[2], 2 );
      real32 const dxMax = -((3*vMax)/(2*dMax[0]))*log( r )*pow((xMax[0]-nodeCoords[a][0])/dMax[0], 2 );
      real32 const dyMax = -((3*vMax)/(2*dMax[1]))*log( r )*pow((xMax[1]-nodeCoords[a][1])/dMax[1], 2 );
      real32 const dzMax = -((3*vMax)/(2*dMax[2]))*log( r )*pow((xMax[2]-nodeCoords[a][2])/dMax[2], 2 );

      if( xMin[0]>nodeCoords[a][0] )
      {
        taperCoeff[a] = LvArray::math::exp( -dxMin*dt );
      }
      else if( nodeCoords[a][0] > xMax[0] )
      {
        taperCoeff[a] = LvArray::math::exp( -dxMax*dt );
      }
      else if( xMin[1]>nodeCoords[a][1] && taperCoeff[a] >= 1.0 )
      {
        taperCoeff[a] = LvArray::math::exp( -dyMin*dt );
      }
      else if( nodeCoords[a][1] > xMax[1] && taperCoeff[a] >= 1.0 )
      {
        taperCoeff[a] = LvArray::math::exp( -dyMax*dt );
      }
      else if( xMin[2]>nodeCoords[a][2] && taperCoeff[a] >= 1.0 )
      {
        taperCoeff[a] = LvArray::math::exp( -dzMin*dt );
      }
      else if( nodeCoords[a][2] > xMax[2] && taperCoeff[a] >= 1.0 )
      {
        taperCoeff[a] = LvArray::math::exp( -dzMax*dt );
      }

    } );
  }


  /**
   * @brief Multiply an array with the taper coefficients
   * @tparam EXEC_POLICY the execution policy
   * @param[in] size the number of nodes
   * @param[in] taperCoeff array which contains the taper coefficient on each node (which will be equal to 1 when we are outside of the
   * taper layers)
   * @param[inout] vector array which is multiplied by the taper array
   */
  template< typename EXEC_POLICY >
  static void
  multiplyByTaperCoeff( localIndex const size,
                        arrayView1d< real32 const > const taperCoeff,
                        arrayView1d< real32 > const vector )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      vector[a] *= taperCoeff[a];
    } );

  }


};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERERNEL_HPP_
