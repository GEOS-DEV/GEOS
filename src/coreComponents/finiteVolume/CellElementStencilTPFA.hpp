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

#ifndef GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "StencilBase.hpp"

namespace geosx
{

/**
 * @struct CellElementStencilTPFA_Traits
 * Struct to predeclare the types and constexpr values of CellElementStencilTPFA so that they may be used in
 * StencilBase.
 */
struct CellElementStencilTPFA_Traits
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
 * @class CellElementStencilTPFAWrapper
 *
 * Class to provide access to the cellElement stencil that may be
 * called from a kernel function.
 */
class CellElementStencilTPFAWrapper : public StencilWrapperBase< CellElementStencilTPFA_Traits >,
  public CellElementStencilTPFA_Traits
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
  CellElementStencilTPFAWrapper( IndexContainerType const & elementRegionIndices,
                                 IndexContainerType const & elementSubRegionIndices,
                                 IndexContainerType const & elementIndices,
                                 WeightContainerType const & weights,
                                 arrayView2d< real64 > const & faceNormal,
                                 arrayView3d< real64 > const & cellToFaceVec,
                                 arrayView1d< real64 > const & transMultiplier )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_faceNormal( faceNormal ),
    m_cellToFaceVec( cellToFaceVec ),
    m_transMultiplier( transMultiplier )
  {}

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
   * Used in ReactiveCompositionalMultiphaseOBL solver for thermal transmissibility computation:
   * here, conductivity is a part of operator and connot be used directly as a coefficient
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weigths w.r.t to the variable
   */
  GEOSX_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
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


private:

  arrayView2d< real64 > m_faceNormal;
  arrayView3d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
};


/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA : public StencilBase< CellElementStencilTPFA_Traits, CellElementStencilTPFA >,
  public CellElementStencilTPFA_Traits
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
                    localIndex const connectorIndex ) override final;

  /**
   * @brief Add the vectors need to compute the transmissiblity to the Stencil.
   * @param[in] transMultiplier the transmissibility multiplier
   * @param[in] faceNormal the normal to the face
   * @param[in] cellToFaceVec distance vector between the cell center and the face
   */
  void addVectors( real64 const & transMultiplier,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[2][3] );

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

  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = CellElementStencilTPFAWrapper;

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

private:
  array2d< real64 > m_faceNormal;
  array3d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;

};

GEOSX_HOST_DEVICE
inline void CellElementStencilTPFAWrapper::computeWeights( localIndex iconn,
                                                           CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                           CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                                                           real64 (& weight)[1][2],
                                                           real64 (& dWeight_dVar )[1][2] ) const
{
  GEOSX_UNUSED_VAR( dCoeff_dVar );

  real64 halfWeight[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i =0; i<2; i++ )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][i];
    localIndex const esr =  m_elementSubRegionIndices[iconn][i];
    localIndex const ei  =  m_elementIndices[iconn][i];

    halfWeight[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3], faceConormal[3];

    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );

    if( LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[er][esr][ei][0], faceNormal );
    halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfWeight[i] < 0.0 )
    {
      LvArray::tensorOps::hadamardProduct< 3 >( faceConormal,
                                                coefficient[er][esr][ei][0],
                                                m_cellToFaceVec[iconn][i] );
      halfWeight[i] = m_weights[iconn][i];
      halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );
    }
  }

  // Do harmonic and arithmetic averaging
  real64 const product = halfWeight[0]*halfWeight[1];
  real64 const sum = halfWeight[0]+halfWeight[1];

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

GEOSX_HOST_DEVICE
inline void CellElementStencilTPFAWrapper::computeWeights( localIndex iconn,
                                                           real64 (& weight)[1][2],
                                                           real64 (& dWeight_dVar )[1][2] ) const
{
  real64 halfWeight[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i =0; i<2; i++ )
  {
    halfWeight[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3];

    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );

    if( LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfWeight[i] < 0.0 )
    {
      halfWeight[i] = m_weights[iconn][i];
      halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], m_cellToFaceVec[iconn][i] );
    }
  }

  // Do harmonic and arithmetic averaging
  real64 const product = halfWeight[0]*halfWeight[1];
  real64 const sum = halfWeight[0]+halfWeight[1];

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


} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
