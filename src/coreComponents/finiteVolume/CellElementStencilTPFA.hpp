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
 * @file CellElementStencilTPFA.hpp
 */

#ifndef GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * Provides access to the cellElement stencil that may be called from a kernel function.
 */
class CellElementStencilTPFAWrapper : public StencilWrapperBase< TwoPointStencilTraits >
{

public:

  //TODO sort out this use
  static constexpr real64 avgWeights = 1.0;


  /**
   * @brief Constructor
   * @param elementRegionIndices The container for the element region indices for each point in each stencil
   * @param elementSubRegionIndices The container for the element sub region indices for each point in each stencil
   * @param elementIndices The container for the element indices for each point in each stencil
   * @param weights The container for the weights for each point in each stencil
   * @param faceNormal Face normal vector
   * @param cellToFaceVec Cell center to face center vector
   * @param transMultiplier Transmissibility multiplier
   * @param geometricStabilizationCoef Geometric coefficient for pressure jump stabilization
   */
  CellElementStencilTPFAWrapper( IndexContainerType const & elementRegionIndices,
                                 IndexContainerType const & elementSubRegionIndices,
                                 IndexContainerType const & elementIndices,
                                 WeightContainerType const & weights,
                                 arrayView2d< real64 > const & faceNormal,
                                 arrayView3d< real64 > const & cellToFaceVec,
                                 arrayView1d< real64 > const & transMultiplier,
                                 arrayView1d< real64 > const & geometricStabilizationCoef );


  /**
   * @brief Compute weights and derivatives w.r.t to one variable based on phase sliced tensor (e.g. diffusion, dispersion)
   * @param[in] iconn connection index
   * @param[in] ip phase index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */

  using StencilWrapperBase< TwoPointStencilTraits >::computeWeights;


  /**
   * @brief Compute weights and derivatives w.r.t to one variable without coefficient
   * Used in ReactiveCompositionalMultiphaseOBL solver for thermal transmissibility computation:
   * here, conductivity is a part of operator and cannot be used directly as a coefficient
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Compute the stabilization weights
   * @param[in] iconn connection index
   * @param[out] stabilizationWeight view weights
   */
  GEOS_HOST_DEVICE
  void computeStabilizationWeights( localIndex iconn,
                                    real64 ( &stabilizationWeight )[1][2] ) const;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  localIndex size() const
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const override
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numPointsInFlux( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

protected:

  GEOS_HOST_DEVICE
  void
  computeWeightsBase( localIndex const iconn,
                      localIndex const (&k)[2],
                      localIndex const icell,
                      arraySlice3d< real64 const > const & coefficient,
                      arraySlice3d< real64 const > const & dCoeff_dVar,
                      real64 & halfWeight,
                      real64 & dHalfWeight_dVar ) const override;

  GEOS_HOST_DEVICE
  void
  computeWeightsBase( localIndex const iconn,
                      localIndex const (&k)[2],
                      localIndex const ielem,
                      arraySlice3d< real64 const >  const & coefficient,
                      arraySlice3d< real64 const >  const & dCoeff_dVar1,
                      arraySlice3d< real64 const >  const & dCoeff_dVar2,
                      real64 & halfWeight,
                      real64 ( & dHalfWeight_dVar )[2] ) const override
  { GEOS_UNUSED_VAR( iconn, k, ielem, coefficient, dCoeff_dVar1, dCoeff_dVar2, halfWeight, dHalfWeight_dVar ); }

  GEOS_HOST_DEVICE
  void
  computeWeightsBase( const geos::localIndex iconn,
                      const geos::localIndex ( & k )[2],
                      const geos::localIndex ielem,
                      const arraySlice3d< const geos::real64 > & coefficient,
                      const arraySlice3d< const geos::real64 > & dCoeff_dVar1,
                      const arraySlice4d< const geos::real64 > & dCoeff_dVar2,
                      real64 & halfWeight,
                      real64 ( & dHalfWeight_dVar ) [2] ) const override
  { GEOS_UNUSED_VAR( iconn, k, ielem, coefficient, dCoeff_dVar1, dCoeff_dVar2, halfWeight, dHalfWeight_dVar ); }

private:

  arrayView2d< real64 > m_faceNormal;
  arrayView3d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
  arrayView1d< real64 > m_geometricStabilizationCoef;


};


/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA final : public StencilBase< TwoPointStencilTraits, CellElementStencilTPFA >
{
public:

  /**
   * @brief Default constructor.
   */
  CellElementStencilTPFA();

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /**
   * @brief Add the vectors need to compute the transmissiblity to the Stencil.
   * @param[in] transMultiplier the transmissibility multiplier
   * @param[in] geometricStabilizationCoef the stabilization weight
   * @param[in] faceNormal the normal to the face
   * @param[in] cellToFaceVec distance vector between the cell center and the face
   */
  void addVectors( real64 const & transMultiplier,
                   real64 const & geometricStabilizationCoef,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[2][3] );

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Reserve the size of the stencil
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size ) override;

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CellElementStencilTPFAWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  array2d< real64 > m_faceNormal;
  array3d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;
  array1d< real64 > m_geometricStabilizationCoef;
};

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::computeWeightsBase( localIndex const iconn,
                                                   localIndex const (&k)[2],
                                                   localIndex const icell,
                                                   arraySlice3d< real64 const > const & coefficient,
                                                   arraySlice3d< real64 const > const & dCoeff_dVar,
                                                   real64 & halfWeight,
                                                   real64 & dHalfWeight_dVar ) const
{
  GEOS_UNUSED_VAR( dCoeff_dVar, dHalfWeight_dVar );

  localIndex const ei = m_elementIndices[iconn][icell];

  halfWeight = m_weights[iconn][icell];

  // Proper computation
  real64 faceNormal[3], cellToFaceVec[3];
  // previously was normalized in container
  LvArray::tensorOps::copy< 3 >( cellToFaceVec, m_cellToFaceVec[iconn][icell] );
  LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

  LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );
  if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
  {
    LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
  }

  real64 faceConormal[3];
  LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[ei][0], faceNormal );
  halfWeight *= LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );

  // correct negative weight issue arising from non-K-orthogonal grids
  if( halfWeight < 0.0 )
  {
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal,
                                              coefficient[ei][0],
                                              cellToFaceVec );
    halfWeight = m_weights[iconn][icell];
    halfWeight *= LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
  }

}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeWeights( localIndex iconn,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{
  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i =0; i<2; i++ )
  {
    halfWeight[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3], cellToFaceVec[3];
    // previously was normalized in container
    LvArray::tensorOps::copy< 3 >( cellToFaceVec, m_cellToFaceVec[iconn][i] );
    LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );

    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfWeight[i] < 0.0 )
    {
      halfWeight[i] = m_weights[iconn][i];
//      halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], m_cellToFaceVec[iconn][i] ); //useless as normalized
// vector it should be 1 always
    }
  }

  averageWeights( iconn, 0 /*connexionIndex*/, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar, 0.5 /*avgCoeff*/ );

}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeStabilizationWeights( localIndex iconn,
                               real64 ( & stabilizationWeight )[1][2] ) const
{
  stabilizationWeight[0][0] = m_geometricStabilizationCoef[iconn];
  stabilizationWeight[0][1] = -m_geometricStabilizationCoef[iconn];
}


} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
