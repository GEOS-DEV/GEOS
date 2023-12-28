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

  /// Coefficient view accessory type
  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::ElementViewConst< VIEWTYPE >;

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

  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       localIndex const ip,
                       CoefficientAccessor< arrayView4d< real64 const > > const &coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const &dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Compute weights and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

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
  localIndex stencilSize( localIndex const index ) const
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

private:

  arrayView2d< real64 > m_faceNormal;
  arrayView3d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
  arrayView1d< real64 > m_geometricStabilizationCoef;

  GEOS_HOST_DEVICE
  void
    averageWeights( const localIndex iconn, real64 ( &weight )[1][2], real64 ( &dWeight_dVar )[1][2],
                    const real64 ( &halfWeight )[2] ) const;

  GEOS_HOST_DEVICE
  void
  computeWeightsBase( localIndex const iconn,
                      localIndex const icell,
                      real64 & halfWeight,
                      arraySlice3d< real64 const > const & coefficient ) const;


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
CellElementStencilTPFAWrapper::
  computeWeights( localIndex const iconn,
                  localIndex const ip,
                  CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{

  GEOS_UNUSED_VAR( dCoeff_dVar );

  real64 halfWeight[2];


  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i = 0; i < 2; ++i )
  {
    localIndex const er = m_elementRegionIndices[iconn][i];
    localIndex const esr = m_elementSubRegionIndices[iconn][i];

    //TODO replace as Sergey LvArray gets merged
    // We are swapping ip phase index and direction to be able to slice properly
    auto coeffNested = coefficient[er][esr];

    LvArray::typeManipulation::CArray< localIndex, 4 > dims, strides;
    dims[0] = coeffNested.dims()[2]; strides[0] = coeffNested.strides()[2];        //swap phase for cell
    dims[1] = coeffNested.dims()[0]; strides[1] = coeffNested.strides()[0];        //increment cell to 2nd pos
    dims[2] = coeffNested.dims()[1]; strides[2] = coeffNested.strides()[1];        //then shift gauss point as well
    dims[3] = coeffNested.dims()[3]; strides[3] = coeffNested.strides()[3];        //direction remain last pos
    ArrayView< real64 const, 4 > coeffSwapped( dims, strides, 0, coeffNested.dataBuffer());

    computeWeightsBase( iconn, i, halfWeight[i], coeffSwapped[ip] );

  }

  // Do harmonic and arithmetic averaging
  averageWeights( iconn, weight, dWeight_dVar, halfWeight );

}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::averageWeights( const localIndex iconn,
                                               real64 (& weight)[1][2],
                                               real64 (& dWeight_dVar)[1][2],
                                               const real64 (& halfWeight)[2] ) const
{
  real64 const product = halfWeight[0] * halfWeight[1];
  real64 const sum = halfWeight[0] + halfWeight[1];

  real64 const harmonicWeight   = sum > 0 ? product / sum : 0.0;
  real64 const arithmeticWeight = sum / 2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member if it is really necessary

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    weight[0][ke] = m_transMultiplier[iconn] * value * (ke == 0 ? 1 : -1);
  }

  dWeight_dVar[0][0] = 0.0;
  dWeight_dVar[0][1] = 0.0;
}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::computeWeightsBase( localIndex const iconn,
                                                   localIndex const icell,
                                                   real64 & halfWeight,
                                                   arraySlice3d< real64 const > const & coefficient ) const
{

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
  computeWeights( localIndex const iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{
  GEOS_UNUSED_VAR( dCoeff_dVar );

  real64 halfWeight[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i = 0; i < 2; ++i )
  {
    localIndex const er = m_elementRegionIndices[iconn][i];
    localIndex const esr = m_elementSubRegionIndices[iconn][i];

    computeWeightsBase( iconn, i, halfWeight[i], coefficient[er][esr] );

  }

  averageWeights( iconn, weight, dWeight_dVar, halfWeight );

}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeWeights( localIndex iconn,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{
  real64 halfWeight[2];

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

  averageWeights( iconn, weight, dWeight_dVar, halfWeight );

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
