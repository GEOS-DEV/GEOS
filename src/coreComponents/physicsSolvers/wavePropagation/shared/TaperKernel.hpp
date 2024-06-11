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
    template< typename EXEC_POLICY>
    static void
    computeTaperCoeff( localIndex const size,
                       arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                       R1Tensor32 const xMin,
                       R1Tensor32 const xMax,
                       R1Tensor32 const dMin,
                       R1Tensor32 const dMax,
                       real32 taperConstant,
                       arrayView1d<real32 > const taperCoeff )
    {
      /// Loop over elements in the subregion, 'l' is the element index within the target set
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
         real32 tmpXmin = (xMin[0]-nodeCoords[a][0])/dMin[0];
         real32 tmpxMax = (nodeCoords[a][0]-xMax[0])/dMax[0];

         real32 tmpYmin = (xMin[1]-nodeCoords[a][1])/dMin[1];
         real32 tmpYmax = (nodeCoords[a][1]-xMax[1])/dMax[1];

         real32 tmpZmin = (xMin[2]-nodeCoords[a][2])/dMin[2];
         real32 tmpZmax = (nodeCoords[a][2]-xMax[2])/dMax[2];

         real32 dist = 0.0;

         if (xMin[0]>nodeCoords[a][0])
         {
           if (xMin[2]>nodeCoords[a][2])
           {
             dist = LvArray::math::min(tmpXmin,tmpZmin);
             if (xMin[1]>nodeCoords[a][1])
             {
                dist= LvArray::math::min(dist,tmpYmin);
             }
             else if (nodeCoords[a][1] > xMax[1])
             {
                dist=LvArray::math::min(dist,tmpYmax);
             }

           }
           else if ( nodeCoords[a][2] > xMax[2])
           {
             dist = LvArray::math::min(tmpXmin,tmpZmax);
             if (xMin[1]>nodeCoords[a][1])
             {
               dist = LvArray::math::min(dist,tmpYmin);
             }
             else if (nodeCoords[a][1] > xMax[1])
             {
                dist=LvArray::math::min(dist,tmpYmax);
             }

           }
           else
           {
            dist=tmpXmin;
            taperCoeff[a] = LvArray::math::exp(-0.015*(dMin[0]-nodeCoords[a][0]));
           }
         }
         else if (nodeCoords[a][0]>xMax[0])
         {
          if (xMin[2]>nodeCoords[a][2])
          {
            dist = LvArray::math::min(tmpxMax,tmpZmin);
            if (xMin[1]>nodeCoords[a][1])
             {
                dist= LvArray::math::min(dist,tmpYmin);
             }
             else if (nodeCoords[a][1] > xMax[1])
             {
                dist=LvArray::math::min(dist,tmpYmax);
             }
          }
          else if ( nodeCoords[a][2] > xMax[2])
           {
             dist = LvArray::math::min(tmpXmin,tmpZmax); //faux
             if (xMin[1]>nodeCoords[a][1])
             {
               dist = LvArray::math::min(dist,tmpYmin);
             }
             else if (nodeCoords[a][1] > xMax[1])
             {
                dist=LvArray::math::min(dist,tmpYmax);
             }

           }
           else
           {
            dist=tmpxMax;
           }
         }
         if (xMin[2]>nodeCoords[a][2])
         {
          if (xMin[1]>nodeCoords[a][1])
          {
             dist= LvArray::math::min(tmpZmin,tmpYmin);
          }
          else if (nodeCoords[a][1] > xMax[1])
          {
             dist=LvArray::math::min(tmpZmin,tmpYmax);
          }
          else
          {
            dist=tmpZmin;
          }

         }
         else if ( nodeCoords[a][2] > xMax[2])
         {
            if (xMin[1]>nodeCoords[a][1])
            {
              dist = LvArray::math::min(tmpZmax,tmpYmin);
            }
            else if (nodeCoords[a][1] > xMax[1])
            {
               dist=LvArray::math::min(tmpZmax,tmpYmax);
            }
            else
            {
              dist=tmpZmax;
            }
         }
         else if (xMin[1]>nodeCoords[a][1])
         {
           dist=tmpYmin;
         }

         else if (nodeCoords[a][1] > xMax[1])
         {
           dist=tmpYmax;
         }
         else
         {
          dist=0.0;

         }

      } );
    }

    template< typename EXEC_POLICY>
    static void
    multiplyByTaperCoeff(localIndex const size,
                         arrayView1d<real32 const > const taperCoeff,
                         arrayView1d<real32 > const vector)
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        vector[a] *= taperCoeff[a];
      } );

    }


  };

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERERNEL_HPP_
