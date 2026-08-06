/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GlobalAdaptationKernels.hpp
 *
 * Residual-based Global Adaptation (GA) indicators for the mixed mimetic discretization.
 * The construction follows four steps:
 *  1) projection of an admissible flow field induced by a nominal uniform gradient,
 *     using the harmonic face-normal diffusive projection across interfaces;
 *  2) evaluation of the localized (cell-wise) TPFA constitutive residual, normalized
 *     by the local pressure drop profile;
 *  3) assembly of the normalized residuals across shared interfaces with respect to
 *     a fixed global orientation of the face normal (Global Adaptation);
 *  4) thresholding to produce the binary stencil activation flag per cell.
 */

#ifndef GEOS_MIXEDMIMETIC_ADAPTIVITY_GLOBALADAPTATIONKERNELS_HPP_
#define GEOS_MIXEDMIMETIC_ADAPTIVITY_GLOBALADAPTATIONKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

namespace mixedMimeticKernels
{

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    case 7:
    { lambda( std::integral_constant< T, 7 >() ); return;}
    case 8:
    { lambda( std::integral_constant< T, 8 >() ); return;}
    case 9:
    { lambda( std::integral_constant< T, 9 >() ); return;}
    case 10:
    { lambda( std::integral_constant< T, 10 >() ); return;}
    case 11:
    { lambda( std::integral_constant< T, 11 >() ); return;}
    case 12:
    { lambda( std::integral_constant< T, 12 >() ); return;}
    case 13:
    { lambda( std::integral_constant< T, 13 >() ); return;}
    default: GEOS_ERROR( GEOS_FMT( "Unknown numFacesInElem value: {}", value ) );
  }
}

} // namespace internal

/******************************** FaceFluxProjectionKernel ********************************/

/**
 * @class FaceFluxProjectionKernel
 * @brief Step 1: on each face, evaluate the discrete normal interface flux induced by a
 *        globally uniform nominal gradient, using the harmonic face-normal diffusive projection.
 */
struct FaceFluxProjectionKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief Launch the face loop computing the projected face fluxes.
   */
  template< typename POLICY >
  static void
  launch( localIndex const numFaces,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const > const & faceCenter,
          arrayView2d< real64 const > const & faceNormal,
          arrayView1d< real64 const > const & faceArea,
          ElementViewConst< arrayView2d< real64 const > > const & elemCenter,
          ElementViewConst< arrayView3d< real64 const > > const & elemPerm,
          real64 const (&gradient)[3],
          real64 const lengthTolerance,
          arrayView1d< real64 > const & projFaceFlux )
  {
    forAll< POLICY >( numFaces, [=] GEOS_HOST_DEVICE ( localIndex const kf )
    {
      real64 sumInvHalfTrans = 0.0;  // sum of d / kappa
      real64 sumDist = 0.0;          // sum of d
      integer elemCounter = 0;

      for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
      {
        localIndex const er  = elemRegionList[kf][k];
        localIndex const esr = elemSubRegionList[kf][k];
        localIndex const ei  = elemList[kf][k];

        bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
        if( onBoundary || !regionFilter.contains( er ) )
        {
          continue;
        }

        // normal half-cell distance from the cell centroid to the face centroid
        real64 cellToFace[3]{};
        for( integer d = 0; d < 3; ++d )
        {
          cellToFace[d] = faceCenter[kf][d] - elemCenter[er][esr][ei][d];
        }
        real64 dist = 0.0;
        for( integer d = 0; d < 3; ++d )
        {
          dist += cellToFace[d] * faceNormal[kf][d];
        }
        dist = LvArray::math::abs( dist );
        dist = LvArray::math::max( dist, lengthTolerance );

        // face-normal directional diffusive component kappa = n . K . n (diagonal K)
        real64 kappa = 0.0;
        real64 kappaScale = 0.0;
        for( integer d = 0; d < 3; ++d )
        {
          kappa += elemPerm[er][esr][ei][0][d] * faceNormal[kf][d] * faceNormal[kf][d];
          kappaScale = LvArray::math::max( kappaScale, elemPerm[er][esr][ei][0][d] );
        }
        // roundoff floor: n.K.n is computed from terms bounded by max_d K_d, so values below
        // eps * max_d K_d are numerically indistinguishable from a collapsed eigenvalue
        kappa = LvArray::math::max( kappa, LvArray::NumericLimits< real64 >::epsilon * kappaScale );

        sumInvHalfTrans += dist / kappa;
        sumDist += dist;
        elemCounter++;
      }

      if( elemCounter == 0 )
      {
        projFaceFlux[kf] = 0.0;
        return;
      }

      // harmonic face-normal diffusive projection (one-sided value on boundary faces)
      real64 const kappaFace = sumDist / sumInvHalfTrans;

      real64 gDotN = 0.0;
      for( integer d = 0; d < 3; ++d )
      {
        gDotN += gradient[d] * faceNormal[kf][d];
      }

      projFaceFlux[kf] = -faceArea[kf] * kappaFace * gDotN;
    } );
  }
};

/******************************** FaceLabelKernel ********************************/

/**
 * @class FaceLabelKernel
 * @brief Classify the faces from the cell-wise stencil activation flags: a face whose
 *        adjacent target cells are all TPFA-compatible gets label 0 (its constitutive row
 *        is exactly diagonal and the flux can be condensed into a two-point expression);
 *        a face adjacent to at least one MFD-compatible cell gets label 1 (live unknown).
 *        When the selected inner product is itself TPFA, the effective operator is diagonal
 *        everywhere and all faces get label 0.
 */
struct FaceLabelKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY >
  static void
  launch( localIndex const numFaces,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          SortedArrayView< localIndex const > const & regionFilter,
          ElementViewConst< arrayView1d< integer const > > const & stencilFlag,
          bool const effectiveTpfa,
          arrayView1d< integer > const & faceStencilLabel )
  {
    forAll< POLICY >( numFaces, [=] GEOS_HOST_DEVICE ( localIndex const kf )
    {
      integer label = 0;
      if( !effectiveTpfa )
      {
        for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[kf][k];
          localIndex const esr = elemSubRegionList[kf][k];
          localIndex const ei  = elemList[kf][k];
          if( er >= 0 && esr >= 0 && ei >= 0 && regionFilter.contains( er ) )
          {
            label = LvArray::math::max( label, stencilFlag[er][esr][ei] );
          }
        }
      }
      faceStencilLabel[kf] = label;
    } );
  }
};

/******************************** LocalResidualKernel ********************************/

/**
 * @class LocalResidualKernel
 * @tparam NF number of faces per element
 * @brief Steps 2-3: evaluate the localized TPFA constitutive residual on the projected pair,
 *        normalize by the local pressure drop profile, and assemble across shared interfaces
 *        with respect to the fixed global orientation of the face normal.
 */
template< integer NF >
struct LocalResidualKernel
{
  template< typename POLICY >
  static void
  launch( localIndex const numElems,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView1d< real64 const > const & elemVolume,
          arrayView3d< real64 const > const & elemPerm,
          arrayView2d< real64 const > const & faceCenter,
          arrayView2d< real64 const > const & faceNormal,
          arrayView1d< real64 const > const & projFaceFlux,
          real64 const (&gradient)[3],
          real64 const lengthTolerance,
          arrayView1d< real64 > const & faceResidual )
  {
    // ||DeltaP_C|| >= |g| * (star-shape radius) on admissible cells: a smaller value signals
    // a degenerate cell, and the normalization saturates at the mesh length tolerance scale
    real64 const dropTolerance = LvArray::tensorOps::l2Norm< 3 >( gradient ) * lengthTolerance;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

      stackArray2d< real64, NF *NF > massMatrix( NF, NF );

      mimeticInnerProduct::TPFAInnerProduct::computeM< NF >( nodePosition,
                                                             faceToNodes,
                                                             elemToFaces[ei],
                                                             elemCenter[ei],
                                                             elemVolume[ei],
                                                             perm,
                                                             lengthTolerance,
                                                             massMatrix );

      real64 localFlux[NF]{};
      real64 pressureDrop[NF]{};
      real64 orientation[NF]{};

      real64 pressureDropNorm = 0.0;
      for( integer i = 0; i < NF; ++i )
      {
        localIndex const kf = elemToFaces[ei][i];

        // sign relating the cell-outward normal to the fixed global face normal
        real64 dotProd = 0.0;
        for( integer d = 0; d < 3; ++d )
        {
          dotProd += ( faceCenter[kf][d] - elemCenter[ei][d] ) * faceNormal[kf][d];
        }
        orientation[i] = ( dotProd >= 0.0 ) ? 1.0 : -1.0;

        // localized (outward-oriented) projected flux
        localFlux[i] = orientation[i] * projFaceFlux[kf];

        // pressure drop profile of the linear probe field: p(x_C) - p(x_f)
        real64 drop = 0.0;
        for( integer d = 0; d < 3; ++d )
        {
          drop += gradient[d] * ( elemCenter[ei][d] - faceCenter[kf][d] );
        }
        pressureDrop[i] = drop;
        pressureDropNorm += drop * drop;
      }
      pressureDropNorm = LvArray::math::sqrt( pressureDropNorm );
      pressureDropNorm = LvArray::math::max( pressureDropNorm, dropTolerance );

      // normalized constitutive residual, assembled on the global face orientation
      for( integer i = 0; i < NF; ++i )
      {
        real64 residual = -pressureDrop[i];
        for( integer j = 0; j < NF; ++j )
        {
          residual += massMatrix( i, j ) * localFlux[j];
        }

        localIndex const kf = elemToFaces[ei][i];
        RAJA::atomicAdd( parallelDeviceAtomic{}, &faceResidual[kf], orientation[i] * residual / pressureDropNorm );
      }
    } );
  }
};

/******************************** MarkingKernel ********************************/

/**
 * @class MarkingKernel
 * @tparam NF number of faces per element
 * @brief Step 4: build the cell-wise indicator as the max of the assembled face residuals
 *        over the faces of the cell, and apply the thresholding criterion.
 */
template< integer NF >
struct MarkingKernel
{
  /**
   * @brief Launch the marking loop.
   * @return the number of locally-owned cells marked as MFD-compatible (ghost cells are
   *         marked too, but excluded from the count so that an MPI sum yields the global count)
   */
  template< typename POLICY >
  static localIndex
  launch( localIndex const numElems,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & faceResidual,
          real64 const tolerance,
          arrayView1d< real64 > const & consistencyIndicator,
          arrayView1d< integer > const & stencilFlag )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, localIndex > numMfdCells( 0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 indicator = 0.0;
      for( integer i = 0; i < NF; ++i )
      {
        indicator = LvArray::math::max( indicator, LvArray::math::abs( faceResidual[elemToFaces[ei][i]] ) );
      }
      consistencyIndicator[ei] = indicator;
      integer const flag = ( indicator > tolerance ) ? 1 : 0;
      stencilFlag[ei] = flag;
      numMfdCells += ( elemGhostRank[ei] < 0 ) ? flag : 0;
    } );

    return numMfdCells.get();
  }
};

} // namespace mixedMimeticKernels

} // namespace geos

#endif //GEOS_MIXEDMIMETIC_ADAPTIVITY_GLOBALADAPTATIONKERNELS_HPP_
