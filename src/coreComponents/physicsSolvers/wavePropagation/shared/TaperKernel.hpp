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
   * @brief Launches the computation of field gradients and divergence for PML region
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] targetSet list of cells in the target set
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemToNodes constant array view of map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] p_n pressure field at time n
   * @param[in] v_n PML auxiliary field at time n
   * @param[in] u_n PML auxiliary field at time n
   * @param[in] xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param[in] xMin PML thickness, left-front-top
   * @param[in] xMax PML thickness, right-back-bottom
   * @param[in] cMin PML wave speed, left-front-top
   * @param[in] cMax PML wave speed, right-back-bottom
   * @param[in] r desired reflectivity of the PML
   * @param[out] grad_n array holding the gradients at time n
   * @param[out] divV_n array holding the divergence at time n
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

      real32 const dx = (3*vMax/2*dMin[0])*log(r)*pow(nodeCoords[a][0]/dMin[0],2);
      real32 const dy = (3*vMax/2*dMin[1])*log(r)*pow(nodeCoords[a][1]/dMin[1],2);
      real32 const dz = (3*vMax/2*dMin[2])*log(r)*pow(nodeCoords[a][2]/dMin[2],2);

      real32 d = dx;

      if(dy < d)
      {
        d = dy;
      }
      if(dz < d)
      {
        d = dz;
      }

      taperCoeff[a] = LvArray::math::exp(-d*dt);

      // if( xMin[0]>nodeCoords[a][0] )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(dMin[0]-nodeCoords[a][0])*(dMin[0]-nodeCoords[a][0]));
      //   taperCoeff[a] = LvArray::math::exp(-dx*dt);
      // }
      // else if( nodeCoords[a][0] > xMax[0] )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(xMax[0]-nodeCoords[a][0])*(xMax[0]-nodeCoords[a][0]));
      // }
      // else if( xMin[1]>nodeCoords[a][1] && taperCoeff[a] >= 1.0 )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(dMin[1]-nodeCoords[a][1])*(dMin[1]-nodeCoords[a][1]));
      // }
      // else if( nodeCoords[a][1] > xMax[1] && taperCoeff[a] >= 1.0 )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(xMax[1]-nodeCoords[a][1])*(xMax[1]-nodeCoords[a][1]));
      // }
      // else if( xMin[2]>nodeCoords[a][2] && taperCoeff[a] >= 1.0 )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(dMin[2]-nodeCoords[a][2])*(dMin[2]-nodeCoords[a][2]));
      // }
      // else if( nodeCoords[a][2] > xMax[2] && taperCoeff[a] >= 1.0 )
      // {
      //   taperCoeff[a] = LvArray::math::exp( -0.0000003*(xMax[2]-nodeCoords[a][2])*(xMax[2]-nodeCoords[a][2]));
      // }

    } );
  }

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
