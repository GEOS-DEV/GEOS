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
 * @file EmbeddedSurfaceToCellStencil.hpp
 */

#ifndef GEOS_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_
#define GEOS_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * @brief Provide access to the EmbeddedSurfaceToCellStencil that may be called from a kernel function.
 */
class EmbeddedSurfaceToCellStencilWrapper : public StencilWrapperBase< TwoPointStencilTraits >
{
public:

  /// Coefficient view accessory type
//  template< typename VIEWTYPE >
//  using CoefficientAccessor = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief Constructor
   * @param elementRegionIndices The container for the element region indices for each point in each stencil
   * @param elementSubRegionIndices The container for the element sub region indices for each point in each stencil
   * @param elementIndices The container for the element indices for each point in each stencil
   * @param weights The container for the weights for each point in each stencil
   */
  EmbeddedSurfaceToCellStencilWrapper( IndexContainerType const & elementRegionIndices,
                                       IndexContainerType const & elementSubRegionIndices,
                                       IndexContainerType const & elementIndices,
                                       WeightContainerType const & weights );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
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
  constexpr localIndex stencilSize( localIndex const index ) const
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
  constexpr localIndex numPointsInFlux( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] ip phase index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
   */

/*  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       localIndex const ip,
                       CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;*/

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
   */
/*  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;*/

    using StencilWrapperBase<TwoPointStencilTraits>::computeWeights;

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
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar1 view accessor to the derivative of the coefficient w.r.t to the variable 1
   * @param[in] dCoeff_dVar2 view accessor to the derivative of the coefficient w.r.t to the variable 2
   * @param[out] weight view weights
   * @param[out] dWeight_dVar1 derivative of the weigths w.r.t to the variable 1
   * @param[out] dWeight_dVar2 derivative of the weigths w.r.t to the variable 2
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar1 )[1][2],
                       real64 ( &dWeight_dVar2 )[1][2] ) const;

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] ip phase index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar1 view accessor to the derivative of the coefficient w.r.t to the variable 1
   * @param[in] dCoeff_dVar2 view accessor to the derivative of the coefficient w.r.t to the variable 2
   * @param[out] weight view weights
   * @param[out] dWeight_dVar1 derivative of the weigths w.r.t to the variable 1
   * @param[out] dWeight_dVar2 derivative of the weigths w.r.t to the variable 2
   */

  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       localIndex const ip,
                       CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar1 )[1][2],
                       real64 ( &dWeight_dVar2 )[1][2] ) const;
  /**
   * @brief Compute the stabilization weights
   * @param[in] iconn connection index
   * @param[out] stabilizationWeight view weights
   */
  GEOS_HOST_DEVICE
  void computeStabilizationWeights( localIndex iconn,
                                    real64 ( & stabilizationWeight )[1][2] ) const
  { GEOS_UNUSED_VAR( iconn, stabilizationWeight ); }

  /**
   * @brief Remove the contribution of the aperture from the weight in the stencil (done before aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void removeHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;

  /**
   * @brief Add the contribution of the aperture to the weight in the stencil (done after aperture update)
   *
   * @param iconn connection index
   * @param hydraulicAperture hydraulic apertures of the fractures
   */
  GEOS_HOST_DEVICE
  void addHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const;

private:

/*  GEOS_HOST_DEVICE
  inline void
    averageWeights2(localIndex const iconn,
                    real64 const ( &halfWeight )[2],
                    real64 const ( &dHalfWeight_dVar )[2][2],
                    real64 ( &weight )[1][2],
                    real64 ( &dWeight_dVar1 )[1][2],
                    real64 ( &dWeight_dVar2 )[1][2] ) const;*/

/*  GEOS_HOST_DEVICE
  inline void
    averageWeights( localIndex const iconn,
                    real64 const ( &halfWeight )[2],
                    real64 const ( &dHalfWeight_dVar )[2],
                    real64 ( &weight )[1][2],
                    real64 ( &dWeight_dVar )[1][2] ) const;*/

  GEOS_HOST_DEVICE
  inline void
    computeWeightsBase( localIndex const iconn,
                        localIndex const ielem,
                        arraySlice3d< real64 const >  const & coefficient,
                        arraySlice3d< real64 const >  const & dCoeff_dVar1,
                        arraySlice3d< real64 const >  const & dCoeff_dVar2,
                        real64 &halfWeight,
                        real64 ( &dHalfWeight_dVar )[2] ) const;

  GEOS_HOST_DEVICE
  inline void
    computeWeightsBase( localIndex const iconn,
                        localIndex const (&k)[2],
                        localIndex const ielem,
                        arraySlice3d< real64 const >  const & coefficient,
                        arraySlice3d< real64 const >  const & dCoeff_dVar,
                        real64 &halfWeight,
                        real64 &dHalfWeight_dVar ) const override;
};

/**
 * @brief Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class EmbeddedSurfaceToCellStencil final : public StencilBase< TwoPointStencilTraits, EmbeddedSurfaceToCellStencil >
{
public:

  virtual void move( LvArray::MemorySpace const space ) override;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = EmbeddedSurfaceToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }



};

/*GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  computeWeights( localIndex iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                  real64 ( & weight )[1][2],
                  real64 ( & dWeight_dVar )[1][2] ) const
{

  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2];

  for( int ielem = 0; ielem < 2; ++ielem )
  {
    localIndex const er = m_elementRegionIndices[iconn][ielem];
    localIndex const esr = m_elementSubRegionIndices[iconn][ielem];
    computeWeightsBase( iconn, ielem, coefficient[er][esr], dCoeff_dVar[er][esr], halfWeight[ielem], dHalfWeight_dVar[ielem] );
  }

    averageWeights(iconn, 0*//*connexionIndex*//*, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar, 1.0 *//*avgCoeff*//* );
}*/

/*GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::computeWeights( const geos::localIndex iconn, const geos::localIndex ip,
                                                     const CoefficientAccessor< arrayView4d< const geos::real64 > > & coefficient,
                                                     const CoefficientAccessor< arrayView4d< const geos::real64 > > & dCoeff_dVar,
                                                     geos::real64 (& weight)[1][2],
                                                     geos::real64 (& dWeight_dVar)[1][2] ) const
{

  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2];


  for( int ielem = 0; ielem < 2; ++ielem )
  {

    // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?
    localIndex const er = m_elementRegionIndices[iconn][ielem];
    localIndex const esr = m_elementSubRegionIndices[iconn][ielem];
    //TODO replace as Sergey LvArray gets merged
    // We are swapping ip phase index and direction to be able to slice properly
    auto coeffNested = coefficient[er][esr];
    auto dCoeffNested_dVar = dCoeff_dVar[er][esr];
    LvArray::typeManipulation::CArray< localIndex, 4 > dims, strides;
    dims[0] = coeffNested.dims()[2];
    strides[0] = coeffNested.strides()[2];                //swap phase for cell
    dims[1] = coeffNested.dims()[0];
    strides[1] = coeffNested.strides()[0];                //increment cell to 2nd pos
    dims[2] = coeffNested.dims()[1];
    strides[2] = coeffNested.strides()[1];                //then shift gauss point as well
    dims[3] = coeffNested.dims()[3];
    strides[3] = coeffNested.strides()[3];                //direction remain last pos
    ArrayView< real64 const, 4 > coeffSwapped( dims, strides, 0, coeffNested.dataBuffer());
    ArrayView< real64 const, 4 > dCoeffSwapped_dVar1( dims, strides, 0, dCoeffNested_dVar.dataBuffer());

    computeWeightsBase( iconn, ielem, coeffSwapped[ip], dCoeffSwapped_dVar1[ip], halfWeight[ielem], dHalfWeight_dVar[ielem] );

  }

    averageWeights(iconn, 0*//*connexionIndex*//*,halfWeight, dHalfWeight_dVar, weight, dWeight_dVar, 1.*//*avgCoeff*//*);
}*/

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  computeWeights( localIndex iconn,
                  real64 ( & weight )[1][2],
                  real64 ( & dWeight_dVar )[1][2] ) const
{
  real64 const t0 = m_weights[iconn][0];
  real64 const t1 = m_weights[iconn][1];

  real64 const sumOfTrans = t0+t1;
  real64 const value = t0*t1/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  dWeight_dVar[0][0] = 0;
  dWeight_dVar[0][1] = 0;
}

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  computeWeights( localIndex const iconn,
                  localIndex const ip,
                  CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                  CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar1 )[1][2],
                  real64 (& dWeight_dVar2) [1][2] ) const
{

  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2][2];


  for( int ielem = 0; ielem < 2; ++ielem )
  {

    // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?
    localIndex const er = m_elementRegionIndices[iconn][ielem];
    localIndex const esr = m_elementSubRegionIndices[iconn][ielem];
    //TODO replace as Sergey LvArray gets merged
    // We are swapping ip phase index and direction to be able to slice properly
    auto coeffNested = coefficient[er][esr];
    auto dCoeffNested_dVar1 = dCoeff_dVar1[er][esr];
    auto dCoeffNested_dVar2 = dCoeff_dVar2[er][esr];
    LvArray::typeManipulation::CArray< localIndex, 4 > dims, strides;
    dims[0] = coeffNested.dims()[2];
    strides[0] = coeffNested.strides()[2];                //swap phase for cell
    dims[1] = coeffNested.dims()[0];
    strides[1] = coeffNested.strides()[0];                //increment cell to 2nd pos
    dims[2] = coeffNested.dims()[1];
    strides[2] = coeffNested.strides()[1];                //then shift gauss point as well
    dims[3] = coeffNested.dims()[3];
    strides[3] = coeffNested.strides()[3];                //direction remain last pos
    ArrayView< real64 const, 4 > coeffSwapped( dims, strides, 0, coeffNested.dataBuffer());
    ArrayView< real64 const, 4 > dCoeffSwapped_dVar1( dims, strides, 0, dCoeffNested_dVar1.dataBuffer());
    ArrayView< real64 const, 4 > dCoeffSwapped_dVar2( dims, strides, 0, dCoeffNested_dVar2.dataBuffer());

    computeWeightsBase( iconn, ielem, coeffSwapped[ip], dCoeffSwapped_dVar1[ip], dCoeffSwapped_dVar2[ip], halfWeight[ielem], dHalfWeight_dVar[ielem] );

  }

    averageWeights(iconn,0/*connexionIndex*/, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar1, dWeight_dVar2, 1.0/*avgCoeff*/);

}



GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  computeWeights( localIndex iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                  CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar1 )[1][2],
                  real64 (& dWeight_dVar2 )[1][2] ) const
{
  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2][2];

  for( int ielem = 0; ielem < 2; ++ielem )
  {
    localIndex const er = m_elementRegionIndices[iconn][ielem];
    localIndex const esr = m_elementSubRegionIndices[iconn][ielem];
    computeWeightsBase( iconn, ielem, coefficient[er][esr], dCoeff_dVar1[er][esr], dCoeff_dVar2[er][esr], halfWeight[ielem], dHalfWeight_dVar[ielem] );
  }

    averageWeights(iconn,/*connexionIndex*/ 0, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar1, dWeight_dVar2, 1.0 /*avgCoeff*/);
}

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::computeWeightsBase( localIndex const iconn,
                                                         localIndex const ielem,
                                                         arraySlice3d< real64 const > const & coefficient,
                                                         arraySlice3d< real64 const > const & dCoeff_dVar1,
                                                         arraySlice3d< real64 const > const & dCoeff_dVar2,
                                                         real64 & halfWeight,
                                                         real64 (& dHalfWeight_dVar)[2] ) const
{


  localIndex const ei = m_elementIndices[iconn][ielem];

  if( ielem == 0 )
  {
    // Will change when implementing collocation points. Will use fracture normal to project the permeability
    halfWeight = m_weights[iconn][ielem] * LvArray::tensorOps::l2Norm< 3 >( coefficient[ei][0] );
    dHalfWeight_dVar[0] = m_weights[iconn][0] * dCoeff_dVar1[ei][0][0];
    dHalfWeight_dVar[1] = m_weights[iconn][0] * dCoeff_dVar2[ei][0][0];
  }
  else
  {
    // We consider the 3rd component of the permeability which is the normal one.
    halfWeight = m_weights[iconn][1] * coefficient[ei][0][2];
    dHalfWeight_dVar[0] = m_weights[iconn][0] * dCoeff_dVar1[ei][0][2];
    dHalfWeight_dVar[1] = m_weights[iconn][0] * dCoeff_dVar2[ei][0][2];
  }

}

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::computeWeightsBase( localIndex const iconn,
                                                         localIndex const (&k)[2],
                                                         localIndex const ielem,
                                                         arraySlice3d< real64 const > const & coefficient,
                                                         arraySlice3d< real64 const > const & dCoeff_dVar,
                                                         real64 & halfWeight,
                                                         real64 &dHalfWeight_dVar ) const
{


  localIndex const ei = m_elementIndices[iconn][ielem];

  if( ielem == 0 )
  {
    // Will change when implementing collocation points. Will use fracture normal to project the permeability
    halfWeight = m_weights[iconn][ielem] * LvArray::tensorOps::l2Norm< 3 >( coefficient[ei][0] );
    dHalfWeight_dVar = m_weights[iconn][ielem] * dCoeff_dVar[ei][0][0];
  }
  else
  {
    // We consider the 3rd component of the permeability which is the normal one.
    halfWeight = m_weights[iconn][ielem] * coefficient[ei][0][2];
    dHalfWeight_dVar = m_weights[iconn][ielem] * dCoeff_dVar[ei][0][2];
  }

}


/*GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::averageWeights2(localIndex const iconn,
                                                     real64 const ( &halfWeight )[2],
                                                     real64 const ( &dHalfWeight_dVar )[2][2],
                                                     real64 ( & weight )[1][2],
                                                     real64 ( & dWeight_dVar1 )[1][2],
                                                     real64 ( & dWeight_dVar2 )[1][2] ) const
{

  //might be used if split differently kept for consistency between wrappers
  GEOS_UNUSED_VAR( iconn );

  real64 const sumOfTrans = halfWeight[0]+halfWeight[1];
  real64 const value = halfWeight[0]*halfWeight[1]/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  // We consider the 3rd component of the permeability which is the normal one.

  dWeight_dVar1[0][0] = ( dHalfWeight_dVar[0][0] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar1[0][1] = ( halfWeight[0] * dHalfWeight_dVar[1][0] * sumOfTrans - dHalfWeight_dVar[1][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

  dWeight_dVar2[0][0] = ( dHalfWeight_dVar[0][1] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar2[0][1] = ( halfWeight[0] * dHalfWeight_dVar[1][1] * sumOfTrans - dHalfWeight_dVar[1][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
}*/

/*
GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::averageWeights( localIndex const iconn,
                                                     real64 const ( &halfWeight )[2],
                                                     real64 const ( &dHalfWeight_dVar )[2],
                                                     real64 ( & weight )[1][2],
                                                     real64 ( & dWeight_dVar )[1][2] ) const
{

  //might be used if split differently kept for consistency between wrappers
  GEOS_UNUSED_VAR( iconn );

  real64 const sumOfTrans = halfWeight[0]+halfWeight[1];
  real64 const value = halfWeight[0]*halfWeight[1]/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  // We consider the 3rd component of the permeability which is the normal one.

  dWeight_dVar[0][0] = ( dHalfWeight_dVar[0] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar[0][1] = ( halfWeight[0] * dHalfWeight_dVar[1] * sumOfTrans - dHalfWeight_dVar[1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

}
*/

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  removeHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er  =  m_elementRegionIndices[iconn][k];
  localIndex const esr =  m_elementSubRegionIndices[iconn][k];
  localIndex const ei  =  m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] * hydraulicAperture[er][esr][ei];
}

GEOS_HOST_DEVICE
inline void
EmbeddedSurfaceToCellStencilWrapper::
  addHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  // only the fracture side is modified, k=1
  localIndex constexpr k = 1;
  localIndex const er  =  m_elementRegionIndices[iconn][k];
  localIndex const esr =  m_elementSubRegionIndices[iconn][k];
  localIndex const ei  =  m_elementIndices[iconn][k];

  m_weights[iconn][k] = m_weights[iconn][k] / hydraulicAperture[er][esr][ei];
}

} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_ */
