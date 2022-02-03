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

#ifndef GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/**
 * @struct EmbeddedSurfaceToCellStencil_Traits
 * Struct to predeclare the types and constexpr values of EmbeddedSurfaceToCellStencil so that they may be used in
 * StencilBase.
 */
struct EmbeddedSurfaceToCellStencil_Traits
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
 * @class EmbeddedSurfaceToCellStencilWrapper
 *
 * Class to provide access to the EmbeddedSurfaceToCellStencil that may be
 * called from a kernel function.
 */
class EmbeddedSurfaceToCellStencilWrapper : public StencilWrapperBase< EmbeddedSurfaceToCellStencil_Traits >,
  public EmbeddedSurfaceToCellStencil_Traits
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
   */
  EmbeddedSurfaceToCellStencilWrapper( IndexContainerType const & elementRegionIndices,
                                       IndexContainerType const & elementSubRegionIndices,
                                       IndexContainerType const & elementIndices,
                                       WeightContainerType const & weights )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights )
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
   * @brief Compute weigths and derivatives w.r.t to one variable without coefficient
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
                       real64 ( &dWeight_dVar2 )[1][2] ) const;


private:

};

/**
 * @class EmbeddedSurfaceToCellStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class EmbeddedSurfaceToCellStencil : public StencilBase< EmbeddedSurfaceToCellStencil_Traits, EmbeddedSurfaceToCellStencil >,
  public EmbeddedSurfaceToCellStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  EmbeddedSurfaceToCellStencil();

  virtual void move( LvArray::MemorySpace const space ) override final;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  /**
   * @brief Add an entry to the stencil.
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] cellCenterToEdgeCenter vectors pointing from the cell center to the edge center
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  void add( localIndex const numPts,
            R1Tensor const * const cellCenterToEdgeCenter,
            localIndex const connectorIndex );


  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = EmbeddedSurfaceToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  StencilWrapper createStencilWrapper() const
  {
    return StencilWrapper( m_elementRegionIndices,
                           m_elementSubRegionIndices,
                           m_elementIndices,
                           m_weights );
  }


  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

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

};

GEOSX_HOST_DEVICE
inline void EmbeddedSurfaceToCellStencilWrapper::computeWeights( localIndex iconn,
                                                                 CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                                 CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                                                                 real64 ( & weight )[1][2],
                                                                 real64 ( & dWeight_dVar )[1][2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * LvArray::tensorOps::l2Norm< 3 >( coefficient[er0][esr0][ei0][0] );
  real64 const t1 = m_weights[iconn][1] * LvArray::tensorOps::l2Norm< 3 >( coefficient[er1][esr1][ei1][0] );

  real64 const sumOfTrans = t0+t1;
  real64 const value = t0*t1/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  real64 const dt0 = m_weights[iconn][0] * dCoeff_dVar[er0][esr0][ei0][0][0];
  real64 const dt1 = m_weights[iconn][1] * dCoeff_dVar[er1][esr1][ei1][0][0];

  dWeight_dVar[0][0] = ( dt0 * t1 * sumOfTrans - dt0 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar[0][1] = ( t0 * dt1 * sumOfTrans - dt1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
}

GEOSX_HOST_DEVICE
inline void EmbeddedSurfaceToCellStencilWrapper::computeWeights( localIndex iconn,
                                                                 real64 ( & weight )[1][2],
                                                                 real64 ( & dWeight_dVar )[1][2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0];
  real64 const t1 = m_weights[iconn][1];

  real64 const sumOfTrans = t0+t1;
  real64 const value = t0*t1/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  real64 const dt0 = m_weights[iconn][0];
  real64 const dt1 = m_weights[iconn][1];

  dWeight_dVar[0][0] = ( dt0 * t1 * sumOfTrans - dt0 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar[0][1] = ( t0 * dt1 * sumOfTrans - dt1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
}

GEOSX_HOST_DEVICE
inline void EmbeddedSurfaceToCellStencilWrapper::computeWeights( localIndex iconn,
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

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * LvArray::tensorOps::l2Norm< 3 >( coefficient[er0][esr0][ei0][0] );
  real64 const t1 = m_weights[iconn][1] * LvArray::tensorOps::l2Norm< 3 >( coefficient[er1][esr1][ei1][0] );

  real64 const sumOfTrans = t0+t1;
  real64 const value = t0*t1/sumOfTrans;

  weight[0][0] = value;
  weight[0][1] = -value;

  real64 const dt0_dVar1 = m_weights[iconn][0] * dCoeff_dVar1[er0][esr0][ei0][0][0];
  real64 const dt1_dVar1 = m_weights[iconn][1] * dCoeff_dVar1[er1][esr1][ei1][0][0];
  real64 const dt0_dVar2 = m_weights[iconn][0] * dCoeff_dVar2[er0][esr0][ei0][0][0];
  real64 const dt1_dVar2 = m_weights[iconn][1] * dCoeff_dVar2[er1][esr1][ei1][0][0];

  dWeight_dVar1[0][0] = ( dt0_dVar1 * t1 * sumOfTrans - dt0_dVar1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar1[0][1] = ( t0 * dt1_dVar1 * sumOfTrans - dt1_dVar1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );

  dWeight_dVar2[0][0] = ( dt0_dVar2 * t1 * sumOfTrans - dt0_dVar2 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar2[0][1] = ( t0 * dt1_dVar2 * sumOfTrans - dt1_dVar2 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
}
} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_ */
