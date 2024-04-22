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
 * @file StencilBase.hpp
 */

#ifndef GEOS_FINITEVOLUME_STENCILBASE_HPP_
#define GEOS_FINITEVOLUME_STENCILBASE_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

/**
 * @brief A collection of properties of a stencil type.
 * @tparam CONTAINER type of container used to store indices and weights
 * @tparam MAX_NUM_POINTS_IN_FLUX maximum number of points connected by a flux
 * @tparam MAX_STENCIL_SIZE maximum number of points in a stencil
 * @tparam MAX_NUM_CONNECTIONS maximum number of connections in a stencil
 */
template< template< typename ... > class CONTAINER,
          localIndex MAX_NUM_POINTS_IN_FLUX,
          localIndex MAX_STENCIL_SIZE,
          localIndex MAX_NUM_CONNECTIONS >
struct StencilTraits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = CONTAINER< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = LvArray::typeManipulation::NestedViewTypeConst< IndexContainerType >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = CONTAINER< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = LvArray::typeManipulation::NestedViewTypeConst< WeightContainerType >;

  /// The array view to type for the stencil weights
  using WeightContainerViewType = LvArray::typeManipulation::NestedViewType< WeightContainerType >;

  /// Maximum number of points the flux
  static constexpr localIndex maxNumPointsInFlux = MAX_NUM_POINTS_IN_FLUX;

  /// Maximum number of points in a stencil
  static constexpr localIndex maxStencilSize = MAX_STENCIL_SIZE;

  /// Maximum number of connections in a stencil
  static constexpr localIndex maxNumConnections = MAX_NUM_CONNECTIONS;
};

/**
 * @brief Describes properties of a standard two-point stencil.
 */
using TwoPointStencilTraits = StencilTraits< array2d, 2, 2, 1 >;

/**
 * @class StencilWrapperBase
 *
 * Class to provide access to the computation of stencil weights that may be
 * called from a kernel function.
 */
template< typename TRAITS >
class StencilWrapperBase : public TRAITS
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
  StencilWrapperBase( typename TRAITS::IndexContainerType const & elementRegionIndices,
                      typename TRAITS::IndexContainerType const & elementSubRegionIndices,
                      typename TRAITS::IndexContainerType const & elementIndices,
                      typename TRAITS::WeightContainerType const & weights ):
    m_elementRegionIndices( elementRegionIndices.toViewConst() ),
    m_elementSubRegionIndices( elementSubRegionIndices.toViewConst() ),
    m_elementIndices( elementIndices.toViewConst() ),
    m_weights( weights.toView() ),
    m_meanPermCoefficient( 1.0 )//as needed in CellElementStencilTPFA and EmbeddedElementStencil
  {};

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOS_HOST_DEVICE
  virtual localIndex stencilSize( localIndex index ) const = 0;

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex index ) const
  {
    return stencilSize( index );
  }
  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementRegionIndices() const { return m_elementRegionIndices; }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementSubRegionIndices() const { return m_elementSubRegionIndices; }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementIndices() const { return m_elementIndices; }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename TRAITS::WeightContainerViewConstType
  getWeights() const { return m_weights; }

  /**
   * @brief Compute weights and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  inline void
    computeWeights( localIndex const iconn,
                    real64 const avgWeight,
                    CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                    CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                    real64 ( &weight )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar )[TRAITS::maxNumConnections][2] ) const;

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
                       real64 const avgWeight,
                       CoefficientAccessor< arrayView4d< real64 const > > const &coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const &dCoeff_dVar,
                       real64 ( &weight )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar )[TRAITS::maxNumConnections][2] ) const;

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
                       real64 const avgWeight,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2] ) const;

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
                       real64 const avgWeight,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2][3] ) const;

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
                       real64 const avgWeight,
                       CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2] ) const;


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
                       real64 const avgWeight,
                       CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                       CoefficientAccessor< arrayView5d< real64 const > > const & dCoeff_dVar2,
                       real64 ( &weight )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                       real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2][3] ) const;

protected:

  /// The container for the element region indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename TRAITS::WeightContainerViewType m_weights;

  /// Mean permeability coefficient
  real64 m_meanPermCoefficient;

  GEOS_HOST_DEVICE
  void
    averageWeights( localIndex const iconn,
                    localIndex const connexionIndex,
                    real64 const ( &halfWeight )[2],
                    real64 const ( &dHalfWeight_dVar )[2],
                    real64 ( &weight )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar )[TRAITS::maxNumConnections][2],
                    real64 const avgCoeff ) const;

  GEOS_HOST_DEVICE
  inline void
    averageWeights( localIndex const iconn,
                    localIndex const connexionIndex,
                    real64 const (&halfWeight)[2],
                    real64 const (&dHalfWeight_dVar)[2][2],
                    real64 ( &weight )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2],
                    real64 const avgCoeff ) const;

  GEOS_HOST_DEVICE
  void
    averageWeights( const localIndex iconn,
                    const localIndex connexionIndex,
                    const real64 ( &halfWeight )[2],
                    const real64 ( &dHalfWeight_dVar )[2][2],
                    real64 ( &weight )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                    real64 ( &dWeight_dVar2 )[TRAITS::maxNumConnections][2][3] ) const;

  GEOS_HOST_DEVICE
  virtual void
  computeWeightsBase( localIndex const iconn,
                      localIndex const (&k)[2],
                      localIndex const icell,
                      arraySlice3d< real64 const > const & coefficient,
                      arraySlice3d< real64 const >  const & dCoeff_dVar,
                      real64 & halfWeight,
                      real64 & dHalfWeight_dVar ) const = 0;

  GEOS_HOST_DEVICE
  virtual void
    computeWeightsBase( localIndex const iconn,
                        localIndex const (&k)[2],
                        localIndex const ielem,
                        arraySlice3d< real64 const >  const & coefficient,
                        arraySlice3d< real64 const >  const & dCoeff_dVar1,
                        arraySlice3d< real64 const >  const & dCoeff_dVar2,
                        real64 &halfWeight,
                        real64 ( &dHalfWeight_dVar )[2] ) const = 0;
  GEOS_HOST_DEVICE
  virtual void
    computeWeightsBase( const geos::localIndex iconn,
                        const geos::localIndex ( &k )[2],
                        const geos::localIndex ielem,
                        const arraySlice3d< const geos::real64 > & coefficient,
                        const arraySlice3d< const geos::real64 > & dCoeff_dVar1,
                        const arraySlice4d< const geos::real64 > & dCoeff_dVar2,
                        real64 & halfWeight,
                        real64 ( &dHalfWeight_dVar ) [2] ) const = 0;



};

GEOS_HOST_DEVICE
template< typename TRAITS >
inline void
StencilWrapperBase< TRAITS >::averageWeights( const geos::localIndex iconn,
                                              const geos::localIndex connexionIndex,
                                              const geos::real64 (& halfWeight)[2],
                                              const geos::real64 (& dHalfWeight_dVar)[2],
                                              geos::real64 (& weight)[TRAITS::maxNumConnections][2],
                                              geos::real64 (& dWeight_dVar)[TRAITS::maxNumConnections][2],
                                              real64 const avgCoeff ) const
{

  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k < 2; ++k )//TODO review hardcoded 2
  {
    //TODO check that
    sumOfTrans += halfWeight[k];
  }

  real64 const arithmeticWeight = avgCoeff * (halfWeight[0]+halfWeight[1]);
  real64 const harmonicWeight   = halfWeight[0]*halfWeight[1] / sumOfTrans;
  real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

  weight[connexionIndex][0] = value;
  weight[connexionIndex][1] = -value;

  real64 dHarmonic[2];
  dHarmonic[0] = ( dHalfWeight_dVar[0] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dHarmonic[1] = ( dHalfWeight_dVar[1] * halfWeight[0] * sumOfTrans - dHalfWeight_dVar[1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

  real64 dArithmetic[2];
  dArithmetic[0] = avgCoeff * dHalfWeight_dVar[0];
  dArithmetic[1] = avgCoeff * dHalfWeight_dVar[1];

  dWeight_dVar[connexionIndex][0] =   ( m_meanPermCoefficient * dHarmonic[0] + (1 - m_meanPermCoefficient) * dArithmetic[0] );
  dWeight_dVar[connexionIndex][1] = -( m_meanPermCoefficient * dHarmonic[1] + (1 - m_meanPermCoefficient) * dArithmetic[1] );

}

GEOS_HOST_DEVICE
template< typename TRAITS >
inline void
StencilWrapperBase< TRAITS >::averageWeights( localIndex const iconn,
                                              localIndex const connexionIndex,
                                              real64 const ( &halfWeight )[2],
                                              real64 const ( &dHalfWeight_dVar )[2][2],
                                              real64 ( & weight )[TRAITS::maxNumConnections][2],
                                              real64 ( & dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                                              real64 ( & dWeight_dVar2 )[TRAITS::maxNumConnections][2],
                                              real64 const avgCoeff ) const
{

  //might be used if split differently kept for consistency between wrappers
  GEOS_UNUSED_VAR( iconn );

  real64 const sumOfTrans = avgCoeff*(halfWeight[0]+halfWeight[1]);
  real64 const value = halfWeight[0]*halfWeight[1]/sumOfTrans;

  weight[connexionIndex][0] = value;
  weight[connexionIndex][1] = -value;

  // We consider the 3rd component of the permeability which is the normal one.

  dWeight_dVar1[connexionIndex][0] = ( dHalfWeight_dVar[0][0] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar1[connexionIndex][1] = ( halfWeight[0] * dHalfWeight_dVar[1][0] * sumOfTrans - dHalfWeight_dVar[1][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

  dWeight_dVar2[connexionIndex][0] = ( dHalfWeight_dVar[0][1] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dWeight_dVar2[connexionIndex][1] = ( halfWeight[0] * dHalfWeight_dVar[1][1] * sumOfTrans - dHalfWeight_dVar[1][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
}
GEOS_HOST_DEVICE
template< typename TRAITS >
inline void
StencilWrapperBase< TRAITS >::averageWeights( const localIndex iconn,
                                              const localIndex connexionIndex,
                                              const real64 (& halfWeight)[2],
                                              const real64 (& dHalfWeight_dVar)[2][2],
                                              real64 (& weight)[TRAITS::maxNumConnections][2],
                                              real64 (& dWeight_dVar1)[TRAITS::maxNumConnections][2],
                                              real64 (& dWeight_dVar2)[TRAITS::maxNumConnections][2][3] ) const
{

  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k<numPointsInFlux( iconn ); ++k )
  {
    //TODO check that
    sumOfTrans += halfWeight[k];
  }

  real64 const arithmeticWeight = 0.25 * (halfWeight[0]+halfWeight[1]);
  real64 const harmonicWeight   = halfWeight[0]*halfWeight[1] / sumOfTrans;
  real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

  weight[connexionIndex][0] = value;
  weight[connexionIndex][1] = -value;

  real64 dHarmonic_dvar1[2];
  dHarmonic_dvar1[0] = ( dHalfWeight_dVar[0][0] * halfWeight[0] * sumOfTrans - dHalfWeight_dVar[0][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dHarmonic_dvar1[1] = ( dHalfWeight_dVar[0][0] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][0] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

  real64 dArithmetic_dvar1[2];
  dArithmetic_dvar1[0] = 0.25 * dHarmonic_dvar1[0];
  dArithmetic_dvar1[1] = 0.25 * dHarmonic_dvar1[1];

  dWeight_dVar1[connexionIndex][0] =    m_meanPermCoefficient * dHarmonic_dvar1[0] + (1 - m_meanPermCoefficient) * dArithmetic_dvar1[0];
  dWeight_dVar1[connexionIndex][1] = -( m_meanPermCoefficient * dHarmonic_dvar1[1] + (1 - m_meanPermCoefficient) * dArithmetic_dvar1[1] );

  real64 dHarmonic_dvar2[2];
  dHarmonic_dvar2[0] = ( dHalfWeight_dVar[0][1] * halfWeight[1] * sumOfTrans - dHalfWeight_dVar[0][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );
  dHarmonic_dvar2[1] = ( halfWeight[0] * dHalfWeight_dVar[1][1] * sumOfTrans - dHalfWeight_dVar[1][1] * halfWeight[0] * halfWeight[1] ) / ( sumOfTrans * sumOfTrans );

  real64 dArithmetic_dvar2[2];
  dArithmetic_dvar2[0] = 0.25 * dHalfWeight_dVar[0][1];
  dArithmetic_dvar2[1] = 0.25 * dHalfWeight_dVar[1][1];

  dWeight_dVar2[connexionIndex][0][0] =   ( m_meanPermCoefficient * dHarmonic_dvar2[0] + (1 - m_meanPermCoefficient) * dArithmetic_dvar2[0] );
  dWeight_dVar2[connexionIndex][1][0] = -( m_meanPermCoefficient * dHarmonic_dvar2[1] + (1 - m_meanPermCoefficient) * dArithmetic_dvar2[1] );

}

GEOS_HOST_DEVICE
template< typename TRAITS >
inline void
StencilWrapperBase< TRAITS >::computeWeights( localIndex const iconn,
                                              real64 const avgWeight,
                                              CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                              CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                                              real64 (& weight)[TRAITS::maxNumConnections][2],
                                              real64 (& dWeight_dVar )[TRAITS::maxNumConnections][2] ) const
{
  GEOS_UNUSED_VAR( dCoeff_dVar );

  localIndex k[2];
  localIndex connectionIndex = 0;

  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1] = k[0] + 1; k[1] < numPointsInFlux( iconn ); ++k[1] )
    {
      //TODO check ordering of k[0] k[1] loops for CellStencil
      real64 halfWeight[2];
      real64 dHalfWeight_dVar[2];

      // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

      for( localIndex ielem = 0; ielem < numPointsInFlux( iconn ); ++ielem )
      {
        localIndex const er = m_elementRegionIndices[iconn][k[ielem]];
        localIndex const esr = m_elementSubRegionIndices[iconn][k[ielem]];

        computeWeightsBase( iconn, k, ielem, coefficient[er][esr], dCoeff_dVar[er][esr], halfWeight[ielem],
                            dHalfWeight_dVar[ielem] );

      }


      averageWeights( iconn, connectionIndex, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar, avgWeight /*avgCoeff*/ );
      connectionIndex++;
    }
  }
}

GEOS_HOST_DEVICE
template< typename TRAITS >
inline void
StencilWrapperBase< TRAITS >::computeWeights( localIndex const iconn,
                                              localIndex const ip,
                                              real64 const avgWeight,
                                              CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                                              CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar,
                                              real64 (& weight)[TRAITS::maxNumConnections][2],
                                              real64 (& dWeight_dVar )[TRAITS::maxNumConnections][2] ) const
{
  localIndex k[2];
  localIndex connectionIndex = 0;

  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1] = k[0] + 1; k[1] < numPointsInFlux( iconn ); ++k[1] )
    {

      real64 halfWeight[2];
      real64 dHalfWeight[2];

      // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

      for( localIndex i = 0; i < 2; ++i )
      {
        localIndex const er = m_elementRegionIndices[iconn][i];
        localIndex const esr = m_elementSubRegionIndices[iconn][i];

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
        ArrayView< real64 const, 4 > coeffSwapped( dims, strides, coeffNested.dataBuffer());
        ArrayView< real64 const, 4 > dCoeffSwapped_dVar( dims, strides, dCoeffNested_dVar.dataBuffer());

        computeWeightsBase( iconn, {0, 1}, i, coeffSwapped[ip], dCoeffSwapped_dVar[ip], halfWeight[i],
                            dHalfWeight[i] );

      }

      // Do harmonic and arithmetic averaging
      averageWeights( iconn, connectionIndex, halfWeight, dHalfWeight, weight, dWeight_dVar,
                      avgWeight );
      connectionIndex++;
    }
  }

}

template< typename TRAITS >
GEOS_HOST_DEVICE
inline void
StencilWrapperBase< TRAITS >::
computeWeights( localIndex const iconn,
                real64 const avgWeight,
                CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                real64 (& weight)[TRAITS::maxNumConnections][2],
                real64 (& dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                real64 (& dWeight_dVar2 )[TRAITS::maxNumConnections][2] ) const
{
  localIndex k[2];
  localIndex connectionIndex = 0;

  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1] = k[0] + 1; k[1] < numPointsInFlux( iconn ); ++k[1] )
    {

      real64 halfWeight[2];
      real64 dHalfWeight_dVar[2][2];

      for( int ielem = 0; ielem < 2; ++ielem )
      {
        localIndex const er = m_elementRegionIndices[iconn][ielem];
        localIndex const esr = m_elementSubRegionIndices[iconn][ielem];
        computeWeightsBase( iconn, ielem, coefficient[er][esr], dCoeff_dVar1[er][esr], dCoeff_dVar2[er][esr],
                            halfWeight[ielem], dHalfWeight_dVar[ielem] );
      }

      averageWeights( iconn, connectionIndex, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar1,
                      dWeight_dVar2, avgWeight );
      connectionIndex++;

    }
  }
}

template< typename TRAITS >
GEOS_HOST_DEVICE
inline void
StencilWrapperBase< TRAITS >::
computeWeights( localIndex const iconn,
                localIndex const ip,
                real64 const avgWeight,
                CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                real64 (& weight)[TRAITS::maxNumConnections][2],
                real64 (& dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                real64 (& dWeight_dVar2) [TRAITS::maxNumConnections][2] ) const
{
  localIndex k[2];
  localIndex connectionIndex = 0;

  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1] = k[0] + 1; k[1] < numPointsInFlux( iconn ); ++k[1] )
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
        strides[0] = coeffNested.strides()[2];                        //swap phase for cell
        dims[1] = coeffNested.dims()[0];
        strides[1] = coeffNested.strides()[0];                        //increment cell to 2nd pos
        dims[2] = coeffNested.dims()[1];
        strides[2] = coeffNested.strides()[1];                        //then shift gauss point as well
        dims[3] = coeffNested.dims()[3];
        strides[3] = coeffNested.strides()[3];                        //direction remain last pos
        ArrayView< real64 const, 4 > coeffSwapped( dims, strides, coeffNested.dataBuffer());
        ArrayView< real64 const, 4 > dCoeffSwapped_dVar1( dims, strides, dCoeffNested_dVar1.dataBuffer());
        ArrayView< real64 const, 4 > dCoeffSwapped_dVar2( dims, strides, dCoeffNested_dVar2.dataBuffer());

        computeWeightsBase( iconn, ielem, coeffSwapped[ip], dCoeffSwapped_dVar1[ip], dCoeffSwapped_dVar2[ip],
                            halfWeight[ielem], dHalfWeight_dVar[ielem] );

      }

      averageWeights( iconn, connectionIndex, halfWeight, dHalfWeight_dVar, weight, dWeight_dVar1, dWeight_dVar2,
                      avgWeight );
      connectionIndex++;
    }
  }

}

template< typename TRAITS >
GEOS_HOST_DEVICE
inline void
StencilWrapperBase< TRAITS >::
computeWeights( localIndex const iconn,
                real64 const avgWeight,
                CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar2,
                real64 ( & weight )[TRAITS::maxNumConnections][2],
                real64 ( & dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                real64 ( & dWeight_dVar2 )[TRAITS::maxNumConnections][2][3] ) const
{
  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {

      real64 halfWeight[2];
      real64 dHalfWeight_dVar[2][2];
      const localIndex elem[2] =  {k[0], k[1]};
      for( localIndex ielem = 0; ielem < 2; ++ielem )
      {

        localIndex const er = m_elementRegionIndices[iconn][k[ielem]];
        localIndex const esr = m_elementSubRegionIndices[iconn][k[ielem]];
        computeWeightsBase( iconn,
                            elem,
                            ielem,
                            coefficient[er][esr],
                            dCoeff_dVar1[er][esr], dCoeff_dVar2[er][esr],
                            halfWeight[ielem], dHalfWeight_dVar[ielem] );

      }

      averageWeights( iconn, connectionIndex,
                      halfWeight, dHalfWeight_dVar,
                      weight, dWeight_dVar1, dWeight_dVar2 );

      connectionIndex++;
    }
  }
}

template< typename TRAITS >
GEOS_HOST_DEVICE
inline void
StencilWrapperBase< TRAITS >::
computeWeights( localIndex const iconn,
                localIndex const ip,
                real64 const avgWeight,
                CoefficientAccessor< arrayView4d< real64 const > > const & coefficient,
                CoefficientAccessor< arrayView4d< real64 const > > const & dCoeff_dVar1,
                CoefficientAccessor< arrayView5d< real64 const > > const & dCoeff_dVar2,
                real64 ( & weight )[TRAITS::maxNumConnections][2],
                real64 ( & dWeight_dVar1 )[TRAITS::maxNumConnections][2],
                real64 ( & dWeight_dVar2 )[TRAITS::maxNumConnections][2][3] ) const
{
  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {

      real64 halfWeight[2];
      real64 dHalfWeight_dVar[2][2];
      const localIndex elem[2] =  {k[0], k[1]};
      for( localIndex ielem = 0; ielem < 2; ++ielem )
      {

        localIndex const er = m_elementRegionIndices[iconn][k[ielem]];
        localIndex const esr = m_elementSubRegionIndices[iconn][k[ielem]];
        //TODO replace as Sergey LvArray gets merged
        // We are swapping ip phase index and direction to be able to slice properly
        auto coeffNested = coefficient[er][esr];
        auto dCoeffNested_dVar1 = dCoeff_dVar1[er][esr];
        auto dCoeffNested_dVar2 = dCoeff_dVar2[er][esr];
        LvArray::typeManipulation::CArray< localIndex, 4 > dims, strides;
        dims[0] = coeffNested.dims()[2];
        strides[0] = coeffNested.strides()[2];                        //swap phase for cell
        dims[1] = coeffNested.dims()[0];
        strides[1] = coeffNested.strides()[0];                        //increment cell to 2nd pos
        dims[2] = coeffNested.dims()[1];
        strides[2] = coeffNested.strides()[1];                        //then shift gauss point as well
        dims[3] = coeffNested.dims()[3];
        strides[3] = coeffNested.strides()[3];                        //direction remain last pos
        ArrayView< real64 const, 4 > coeffSwapped( dims, strides, coeffNested.dataBuffer());
        ArrayView< real64 const, 4 > dCoeffSwapped_dVar1( dims, strides, dCoeffNested_dVar1.dataBuffer());

        LvArray::typeManipulation::CArray< localIndex, 5 > ddims, dstrides;
        ddims[0] = dCoeffNested_dVar2.dims()[2];
        dstrides[0] = dCoeffNested_dVar2.strides()[2];                        //swap phase for cell
        ddims[1] = dCoeffNested_dVar2.dims()[0];
        dstrides[1] = dCoeffNested_dVar2.strides()[0];                        //increment cell to 2nd pos
        ddims[2] = dCoeffNested_dVar2.dims()[1];
        dstrides[2] = dCoeffNested_dVar2.strides()[1];                        //then shift gauss point as well
        ddims[3] = dCoeffNested_dVar2.dims()[3];
        dstrides[3] = dCoeffNested_dVar2.strides()[3];
        ddims[4] = dCoeffNested_dVar2.dims()[4];
        dstrides[4] = dCoeffNested_dVar2.strides()[4];
        ArrayView< real64 const, 5 > dCoeffSwapped_dVar2( ddims, dstrides, dCoeffNested_dVar2.dataBuffer());

        computeWeightsBase( iconn,
                            elem,
                            ielem,
                            coefficient[er][esr],
                            dCoeff_dVar1[er][esr], dCoeff_dVar2[er][esr],
                            halfWeight[ielem], dHalfWeight_dVar[ielem] );

      }

      averageWeights( iconn, connectionIndex,
                      halfWeight, dHalfWeight_dVar,
                      weight, dWeight_dVar1, dWeight_dVar2 );

      connectionIndex++;
    }
  }
}


/**
 * @brief Provides management of the interior stencil points when using Two-Point flux approximation.
 * @tparam TRAITS the traits class describing properties of the stencil
 * @tparam LEAFCLASS derived type for CRTP
 */
template< typename TRAITS, typename LEAFCLASS >
class StencilBase : public TRAITS
{
public:

  /**
   * @brief Destructor.
   */
  virtual ~StencilBase() = default;

  /**
   * @brief Reserve the size of the stencil.
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size );

  /**
   * @brief Move the data arrays associated with the stencil to a specified
   *   memory space.
   * @param space The target memory space.
   *
   * @note The existence of this function indicates we need to redesign the
   * stencil classes.
   */
  virtual void move( LvArray::MemorySpace const space );


  /**
   * @brief Add an entry to the stencil.
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] elementRegionIndices The element region indices for each point in the stencil entry
   * @param[in] elementSubRegionIndices The element sub-region indices for each point in the stencil entry
   * @param[in] elementIndices The element indices for each point in the stencil entry
   * @param[in] weights The weights each point in the stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) = 0;

  /**
   * @brief Zero weights for a stencil entry.
   * @param[in] connectorIndex The index of the connector element that the stencil acts across for which the weights are
   *                           to be zero.
   * @return True if a valid connectorIndex was found, and had its corresponding weights set to zero.
   */
  virtual bool zero( localIndex const connectorIndex );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const = 0;

  /**
   * @brief Set the name used in data movement logging callbacks.
   * @param name the name prefix for the stencil's data arrays
   */
  void setName( string const & name );

  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementRegionIndices() const { return m_elementRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementSubRegionIndices() const { return m_elementSubRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementIndices() const { return m_elementIndices.toViewConst(); }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename TRAITS::WeightContainerViewConstType
  getWeights() const { return m_weights.toViewConst(); }

protected:

  /// The container for the element region indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename TRAITS::WeightContainerType m_weights;

  /// The map that provides the stencil index given the index of the underlying connector object.
  unordered_map< localIndex, localIndex > m_connectorIndices;
};



template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 2 );
  m_elementSubRegionIndices.reserve( size * 2 );
  m_elementIndices.reserve( size * 2 );
  m_weights.reserve( size * 2 );
}


template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
bool StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::zero( localIndex const connectorIndex )
{
  return
    executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for( localIndex i = 0; i < static_cast< LEAFCLASS * >(this)->stencilSize( connectorIndex ); ++i )
    {
      m_weights[connectionListIndex][i] = 0;
    }
  } );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::setName( string const & name )
{
  m_elementRegionIndices.setName( name + "/elementRegionIndices" );
  m_elementSubRegionIndices.setName( name + "/elementSubRegionIndices" );
  m_elementIndices.setName( name + "/elementIndices" );
  m_weights.setName( name + "/weights" );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::move( LvArray::MemorySpace const space )
{
  m_elementRegionIndices.move( space, true );
  m_elementSubRegionIndices.move( space, true );
  m_elementIndices.move( space, true );
  m_weights.move( space, true );
}

} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_STENCILBASE_HPP_ */
