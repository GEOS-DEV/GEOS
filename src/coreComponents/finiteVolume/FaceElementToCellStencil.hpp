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
 * @file FaceElementToCellStencil.hpp
 */

#ifndef GEOS_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_
#define GEOS_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * @brief Provides access to the FaceElementToCellStencil that may be called from a kernel function.
 */
class FaceElementToCellStencilWrapper : public StencilWrapperBase< TwoPointStencilTraits >
{
public:

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
   */
  FaceElementToCellStencilWrapper( IndexContainerType const & elementRegionIndices,
                                   IndexContainerType const & elementSubRegionIndices,
                                   IndexContainerType const & elementIndices,
                                   WeightContainerType const & weights,
                                   arrayView2d< real64 > const & faceNormal,
                                   arrayView2d< real64 > const & cellToFaceVec,
                                   arrayView1d< real64 > const & transMultiplier );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex size() const
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const override
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  using StencilWrapperBase< TwoPointStencilTraits >::computeWeights;

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable without coefficient
   * Used in ReactiveCompositionalMultiphaseOBL solver for thermal transmissibility computation:
   * here, conductivity is a part of operator and connot be used directly as a coefficient
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
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
                                    real64 ( & stabilizationWeight )[1][2] ) const
  {
    GEOS_UNUSED_VAR( iconn, stabilizationWeight );
  }

  /**
   * @brief Remove the contribution of the aperture from the weight in the stencil (done before aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void removeHydraulicApertureContribution( localIndex const iconn,
                                            ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;

  /**
   * @brief Add the contribution of the aperture to the weight in the stencil (done after aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void addHydraulicApertureContribution( localIndex const iconn,
                                         ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;

private:

  /// Face normal vector
  arrayView2d< real64 > m_faceNormal;

  /// Cell center to face center vector
  arrayView2d< real64 > m_cellToFaceVec;

  /// Transmissibility multiplier
  arrayView1d< real64 > m_transMultiplier;


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
                        localIndex const icell,
                        arraySlice3d< real64 const > const & coefficient,
                        arraySlice3d< real64 const > const & dCoeff_dVar1,
                        arraySlice3d< real64 const > const & dCoeff_dVar2,
                        real64 & halfWeight,
                        real64 ( &dHalfWeight_dVar )[2] ) const override;


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

};

/**
 * @class FaceElementToCellStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class FaceElementToCellStencil final : public StencilBase< TwoPointStencilTraits, FaceElementToCellStencil >
{
public:

  /**
   * @brief Default constructor.
   */
  FaceElementToCellStencil();

  virtual void move( LvArray::MemorySpace const space ) override;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /**
   * @brief Adds the vectors need to compute weights needed in kernels
   * @param transMultiplier transmissibility multiplier
   * @param faceNormal face normal vector
   * @param cellToFaceVec cell to face vector
   */
  void addVectors( real64 const & transMultiplier,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[3] );

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FaceElementToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override { return m_elementRegionIndices.size( 0 ); }

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

private:

  array2d< real64 > m_faceNormal;
  array2d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;
};


GEOS_HOST_DEVICE
inline void
FaceElementToCellStencilWrapper::computeWeightsBase( localIndex const iconn,
                                                     localIndex const (&k)[2],
                                                     localIndex const ielem,
                                                     arraySlice3d< real64 const > const & coefficient,
                                                     arraySlice3d< real64 const > const & dCoeff_dVar,
                                                     real64 & halfWeight,
                                                     real64 & dHalfWeight_dVar ) const
{

  GEOS_UNUSED_VAR( k );

  localIndex const ei = m_elementIndices[iconn][ielem];


  if( ielem == 0 )
  {
    real64 faceConormal[3];
    // Will change when implementing collocation points. Will use fracture normal to project the permeability
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[ei][0], m_faceNormal[iconn] );
    halfWeight = m_weights[iconn][ielem] *
                 LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal );
    dHalfWeight_dVar = m_weights[iconn][ielem] * dCoeff_dVar[ei][0][0];
  }
  else        //ielem == 1
  {
    // We consider the 3rd component of the permeability which is the normal one.
    halfWeight = m_weights[iconn][ielem] * coefficient[ei][0][2];
    dHalfWeight_dVar = m_weights[iconn][ielem] * dCoeff_dVar[ei][0][2];
  }

}

GEOS_HOST_DEVICE
inline void
FaceElementToCellStencilWrapper::computeWeightsBase( localIndex const iconn,
                                                     localIndex const (&k)[2],
                                                     localIndex const ielem,
                                                     arraySlice3d< real64 const > const & coefficient,
                                                     arraySlice3d< real64 const > const & dCoeff_dVar1,
                                                     arraySlice3d< real64 const > const & dCoeff_dVar2,
                                                     real64 & halfWeight,
                                                     real64 (& dHalfWeight_dVar)[2] ) const
{
  GEOS_UNUSED_VAR( k );
  localIndex const ei = m_elementIndices[iconn][ielem];
  if( ielem == 0 )
  {
    real64 faceConormal[3];
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[ei][0], m_faceNormal[iconn] );
    // Will change when implementing collocation points. Will use fracture normal to project the permeability
    halfWeight =m_transMultiplier[iconn] * m_weights[iconn][0] * LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal );
    dHalfWeight_dVar[0] = m_transMultiplier[iconn] * m_weights[iconn][0] * dCoeff_dVar1[ei][0][0];
    dHalfWeight_dVar[1] = m_transMultiplier[iconn] * m_weights[iconn][0] * dCoeff_dVar2[ei][0][0];
  }
  else
  {
    // We consider the 3rd component of the permeability which is the normal one.
    halfWeight = m_transMultiplier[iconn] * m_weights[iconn][1] * coefficient[ei][0][2];
    dHalfWeight_dVar[0] = m_transMultiplier[iconn] * m_weights[iconn][0] * dCoeff_dVar1[ei][0][2];
    dHalfWeight_dVar[1] = m_transMultiplier[iconn] * m_weights[iconn][0] * dCoeff_dVar2[ei][0][2];
  }

}

GEOS_HOST_DEVICE
inline void
FaceElementToCellStencilWrapper::computeWeights( localIndex iconn,
                                                 real64 ( & weight )[1][2],
                                                 real64 ( & dWeight_dVar )[1][2] ) const
{
  // Will change when implementing collocation points.
  real64 const t0 =
    m_weights[iconn][0] * LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], m_faceNormal[iconn] );
  real64 const t1 = m_weights[iconn][1];

  real64 const sumOfTrans = t0 + t1;
  real64 const value = m_transMultiplier[iconn] * t0 * t1 / sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  dWeight_dVar[0][0] = 0.0;
  dWeight_dVar[0][1] = 0.0;
}

GEOS_HOST_DEVICE
inline void
FaceElementToCellStencilWrapper::
  removeHydraulicApertureContribution( localIndex const iconn,
                                       ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er = m_elementRegionIndices[iconn][k];
  localIndex const esr = m_elementSubRegionIndices[iconn][k];
  localIndex const ei = m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] * hydraulicAperture[er][esr][ei];
}

GEOS_HOST_DEVICE
inline void
FaceElementToCellStencilWrapper::
  addHydraulicApertureContribution( localIndex const iconn,
                                    ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er = m_elementRegionIndices[iconn][k];
  localIndex const esr = m_elementSubRegionIndices[iconn][k];
  localIndex const ei = m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] / hydraulicAperture[er][esr][ei];
}

} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_ */
