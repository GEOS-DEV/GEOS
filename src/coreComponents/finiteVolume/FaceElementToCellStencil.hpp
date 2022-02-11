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

#ifndef GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/**
 * @struct FaceElementToCellStencil_Traits
 * Struct to predeclare the types and constexpr values of FaceElementToCellStencil so that they may be used in
 * StencilBase.
 */
struct FaceElementToCellStencil_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = array2d< localIndex >;

  /// The array view type for the stencil indices
  using IndexContainerViewType = arrayView2d< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = arrayView2d< localIndex const >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = array2d< real64 >;

  /// The array view type for the stencil weights
  using WeightContainerViewType = arrayView2d< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = arrayView2d< real64 const >;

  /// Number of points the flux is between (always 2 for TPFA)
  static constexpr localIndex NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil (this is 2 for TPFA)
  static constexpr localIndex MAX_STENCIL_SIZE = 2;

  /// Maximum number of connections in a stencil
  static constexpr localIndex MAX_NUM_OF_CONNECTIONS = 1;
};

/**
 * @class FaceElementToCellStencilWrapper
 *
 * Class to provide access to the FaceElementToCellStencil that may be
 * called from a kernel function.
 */
class FaceElementToCellStencilWrapper : public StencilWrapperBase< FaceElementToCellStencil_Traits >,
  public FaceElementToCellStencil_Traits
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
   */
  FaceElementToCellStencilWrapper( IndexContainerType const & elementRegionIndices,
                                   IndexContainerType const & elementSubRegionIndices,
                                   IndexContainerType const & elementIndices,
                                   WeightContainerType const & weights,
                                   arrayView2d< real64 > const & faceNormal,
                                   arrayView2d< real64 > const & cellToFaceVec,
                                   arrayView1d< real64 > const & transMultiplier )
    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_faceNormal( faceNormal ),
    m_cellToFaceVec( cellToFaceVec ),
    m_transMultiplier( transMultiplier )
  {}

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return MAX_STENCIL_SIZE;
  }

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex numPointsInFlux( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return NUM_POINT_IN_FLUX;
  }

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
   */
  GEOSX_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
   */
  GEOSX_HOST_DEVICE
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
  GEOSX_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  dCoeff_dVar1,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  dCoeff_dVar2,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar1 )[1][2],
                       real64 ( &dWeight_dVar2 )[1] [2] ) const;

private:

  /// Face normal vector
  arrayView2d< real64 > m_faceNormal;

  /// Cell center to face center vector
  arrayView2d< real64 > m_cellToFaceVec;

  /// Transmissibility multiplier
  arrayView1d< real64 > m_transMultiplier;
};

/**
 * @class FaceElementToCellStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class FaceElementToCellStencil : public StencilBase< FaceElementToCellStencil_Traits, FaceElementToCellStencil >,
  public FaceElementToCellStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  FaceElementToCellStencil();

  virtual void move( LvArray::MemorySpace const space ) override final;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

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
  using StencilWrapper = FaceElementToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  StencilWrapper createStencilWrapper() const
  {
    return StencilWrapper( m_elementRegionIndices,
                           m_elementSubRegionIndices,
                           m_elementIndices,
                           m_weights,
                           m_faceNormal,
                           m_cellToFaceVec,
                           m_transMultiplier );
  }


  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Reserve the size of the stencil
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size ) override final;

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return MAX_STENCIL_SIZE;
  }

private:
  array2d< real64 > m_faceNormal;
  array2d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;

};

GEOSX_HOST_DEVICE
inline void FaceElementToCellStencilWrapper::computeWeights( localIndex iconn,
                                                             CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                             CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                                                             real64 ( & weight )[1][2],
                                                             real64 ( & dWeight_dVar )[1][2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  real64 halfWeight = m_weights[iconn][0];

  real64 faceConormal[3];

  LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[er0][esr0][ei0][0], m_faceNormal[iconn] );
  halfWeight *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal ) * m_transMultiplier[iconn];

  weight[0][0] = halfWeight;
  weight[0][1] = -halfWeight;

  dWeight_dVar[0][0] = 0.0 * dCoeff_dVar[er0][esr0][ei0][0][0];
  dWeight_dVar[0][1] = 0.0;
}

GEOSX_HOST_DEVICE
inline void FaceElementToCellStencilWrapper::computeWeights( localIndex iconn,
                                                             real64 ( & weight )[1][2],
                                                             real64 ( & dWeight_dVar )[1][2] ) const
{
  real64 halfWeight = m_weights[iconn][0];

  halfWeight *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], m_faceNormal[iconn] ) * m_transMultiplier[iconn];

  weight[0][0] = halfWeight;
  weight[0][1] = -halfWeight;

  dWeight_dVar[0][0] = 0.0;
  dWeight_dVar[0][1] = 0.0;
}

GEOSX_HOST_DEVICE
inline void FaceElementToCellStencilWrapper::computeWeights( localIndex iconn,
                                                             CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                             CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                                                             CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                                                             real64 (& weight)[1][2],
                                                             real64 (& dWeight_dVar1 )[1][2],
                                                             real64 (& dWeight_dVar2 )[1][2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  real64 halfWeight = m_weights[iconn][0];

  real64 faceConormal[3];

  LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[er0][esr0][ei0][0], m_faceNormal[iconn] );
  halfWeight *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal ) * m_transMultiplier[iconn];

  weight[0][0] = halfWeight;
  weight[0][1] = -halfWeight;

  dWeight_dVar1[0][0] = 0.0 * dCoeff_dVar1[er0][esr0][ei0][0][0];
  dWeight_dVar1[0][1] = 0.0;

  dWeight_dVar2[0][0] = 0.0 * dCoeff_dVar2[er0][esr0][ei0][0][0];
  dWeight_dVar2[0][1] = 0.0;
}


} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_ */
