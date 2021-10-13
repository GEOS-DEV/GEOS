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
 * @file SurfaceElementStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{

/// @cond DO_NOT_DOCUMENT
// TODO remove! This option allows for the creation of new mass inside a newly
// created FaceElement. The new mass will be equal to:
// creationMass = defaultDensity * defaultAperture * faceArea.
// If 0, then the beginning of step density is artificially set to zero...which
// may cause some newton convergence problems.
#define ALLOW_CREATION_MASS 1


// TODO remove! This option sets the pressure in a newly created FaceElement to
// be the lowest value of all attached non-new FaceElements.
#define SET_CREATION_PRESSURE 1

// TODO remove! This option sets the nodal displacements attached a newly
// created FaceElement to some scalar fraction of the aperture of the
// lowest attached non-new FaceElements.
#define SET_CREATION_DISPLACEMENT 0
/// @endcond

/**
 * @struct SurfaceElementStencil_Traits
 * Struct to predeclare the types and constexpr values of SurfaceElementStencil so that they may be used in
 * StencilBase.
 */
struct SurfaceElementStencil_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = ArrayOfArrays< localIndex >;

  /// The array view type for the stencil indices
  using IndexContainerViewType = ArrayOfArraysView< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = ArrayOfArraysView< localIndex const >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = ArrayOfArrays< real64 >;

  /// The array view type for the stencil weights
  using WeightContainerViewType = ArrayOfArraysView< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = ArrayOfArraysView< real64 const >;

  /// Number of points the flux is between (normally 2)
  static localIndex constexpr NUM_POINT_IN_FLUX = 6;

  /// Maximum number of points in a stencil
  static localIndex constexpr MAX_STENCIL_SIZE = 6;

  /// Maximum number of connections in a stencil
  static localIndex constexpr MAX_NUM_OF_CONNECTIONS = MAX_STENCIL_SIZE * (MAX_STENCIL_SIZE - 1) / 2;
};

/**
 * @class SurfaceElementStencilWrapper
 *
 * Class to provide access to the SurfaceElementStencil that may be
 * called from a kernel function.
 */
class SurfaceElementStencilWrapper : public StencilWrapperBase< SurfaceElementStencil_Traits >,
  public SurfaceElementStencil_Traits
{
public:

  /// Threshold for the application of the permeability multiplier
  static constexpr real64 MULTIPLIER_THRESHOLD = 1e-10;

  /// Coefficient view accessory type
  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief Constructor
   * @param elementRegionIndices The container for the element region indices for each point in each stencil
   * @param elementSubRegionIndices The container for the element sub region indices for each point in each stencil
   * @param elementIndices The container for the element indices for each point in each stencil
   * @param weights The container for the weights for each point in each stencil
   * @param cellCenterToEdgeCenters Cell center to Edge center vector
   * @param meanPermCoefficient Mean permeability coefficient
   */
  SurfaceElementStencilWrapper( IndexContainerType const & elementRegionIndices,
                                IndexContainerType const & elementSubRegionIndices,
                                IndexContainerType const & elementIndices,
                                WeightContainerType const & weights,
                                ArrayOfArrays< R1Tensor > const & cellCenterToEdgeCenters,
                                real64 const meanPermCoefficient )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_cellCenterToEdgeCenters( cellCenterToEdgeCenters.toView() ),
    m_meanPermCoefficient( meanPermCoefficient )
  {}

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }


  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex numPointsInFlux( localIndex index ) const
  {
    return stencilSize( index );
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
                       real64 ( &weight )[MAX_NUM_OF_CONNECTIONS][2],
                       real64 ( &dWeight_dVar )[MAX_NUM_OF_CONNECTIONS][2] ) const;


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
                       real64 ( &weight )[MAX_NUM_OF_CONNECTIONS][2],
                       real64 ( &dWeight_dVar1 )[MAX_NUM_OF_CONNECTIONS][2],
                       real64 ( &dWeight_dVar2 )[MAX_NUM_OF_CONNECTIONS][2] ) const;

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] coefficientMultiplier view accessor to the coefficient multiplier used to compute the weights
   * @param[in] gravityVector gravity vector
   * @param[out] weight view weights
   */
  GEOSX_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficientMultiplier,
                       R1Tensor const & gravityVector,
                       real64 ( &weight )[MAX_NUM_OF_CONNECTIONS][2] ) const;

  /**
   * @brief Compute weigths and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient1 view accessor to the first coefficient used to compute the first weights
   * @param[in] coefficient1Multiplier view accessor to the coefficient multiplier used to compute the first weights
   * @param[in] coefficient2 view accessor to the first coefficient used to compute the second weights
   * @param[in] gravityVector gravity vector
   * @param[out] weight1 view on the first weights
   * @param[out] weight2 view on the second weights
   * @param[out] geometricWeight view on the purely geometric weights
   */
  GEOSX_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficient1,
                       CoefficientAccessor< arrayView3d< real64 const > > const &  coefficient1Multiplier,
                       CoefficientAccessor< arrayView1d< real64 const > > const &  coefficient2,
                       R1Tensor const & gravityVector,
                       real64 ( &weight1 )[NUM_POINT_IN_FLUX],
                       real64 ( &weight2 )[NUM_POINT_IN_FLUX],
                       real64 ( &geometricWeight )[NUM_POINT_IN_FLUX] ) const;


  /**
   * @brief Accessor to the CellCenterToEdgeCenter vector
   * @return the view const to the CellCenterToEdgeCenter vector
   */
  ArrayOfArraysView< R1Tensor const > getCellCenterToEdgeCenters() const
  { return m_cellCenterToEdgeCenters.toViewConst(); }

private:
  /// Cell center to Edge center vector
  ArrayOfArraysView< R1Tensor > m_cellCenterToEdgeCenters;

  /// Mean permeability coefficient
  real64 m_meanPermCoefficient;
};

/**
 * @class SurfaceElementStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class SurfaceElementStencil : public StencilBase< SurfaceElementStencil_Traits, SurfaceElementStencil >,
  public SurfaceElementStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  SurfaceElementStencil();

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
  using StencilWrapper = SurfaceElementStencilWrapper;

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
                           m_cellCenterToEdgeCenters,
                           m_meanPermCoefficient );
  }


  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }

  /**
   * @brief Give the array of vectors pointing from the cell center to the edge center.
   * @return The array of vectors pointing from the cell center to the edge center
   */
  ArrayOfArraysView< R1Tensor const > getCellCenterToEdgeCenters() const
  { return m_cellCenterToEdgeCenters.toViewConst(); }

  /**
   * @brief sets the value of the mean perm conefficient
   * @param meanPermCoefficient value to be set
   */
  void setMeanPermCoefficient( real64 const & meanPermCoefficient )
  {
    m_meanPermCoefficient = meanPermCoefficient;
  }

private:

  /// Distance between the center of the face element and the center of the connecting edge.
  ArrayOfArrays< R1Tensor > m_cellCenterToEdgeCenters;

  /// Mean permeability coefficient
  real64 m_meanPermCoefficient;

};

GEOSX_HOST_DEVICE
inline void SurfaceElementStencilWrapper::computeWeights( localIndex iconn,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                                                          real64 ( & weight )[MAX_NUM_OF_CONNECTIONS][2],
                                                          real64 ( & dWeight_dVar )[MAX_NUM_OF_CONNECTIONS][2] ) const
{

  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k<numPointsInFlux( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    sumOfTrans += coefficient[er][esr][ei][0][0] * m_weights[iconn][k];
  }

  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {
      localIndex const er0  =  m_elementRegionIndices[iconn][k[0]];
      localIndex const esr0 =  m_elementSubRegionIndices[iconn][k[0]];
      localIndex const ei0  =  m_elementIndices[iconn][k[0]];

      localIndex const er1  =  m_elementRegionIndices[iconn][k[1]];
      localIndex const esr1 =  m_elementSubRegionIndices[iconn][k[1]];
      localIndex const ei1  =  m_elementIndices[iconn][k[1]];

      real64 const t0 = m_weights[iconn][0] * coefficient[er0][esr0][ei0][0][0]; // this is a bit insane to access perm
      real64 const t1 = m_weights[iconn][1] * coefficient[er1][esr1][ei1][0][0];

      real64 const harmonicWeight   = t0*t1 / sumOfTrans;
      real64 const arithmeticWeight = 0.25 * (t0+t1);

      real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

      weight[connectionIndex][0] = value;
      weight[connectionIndex][1] = -value;

      real64 const dt0 = m_weights[iconn][0] * dCoeff_dVar[er0][esr0][ei0][0][0];
      real64 const dt1 = m_weights[iconn][1] * dCoeff_dVar[er1][esr1][ei1][0][0];

      real64 dHarmonic[2];
      dHarmonic[0] = ( dt0 * t1 * sumOfTrans - dt0 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
      dHarmonic[1] = ( t0 * dt1 * sumOfTrans - dt1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );

      real64 dArithmetic[2];
      dArithmetic[0] = 0.25 * dt0;
      dArithmetic[1] = 0.25 * dt1;

      dWeight_dVar[connectionIndex][0] = m_meanPermCoefficient * dHarmonic[0] + (1 - m_meanPermCoefficient) * dArithmetic[0];
      dWeight_dVar[connectionIndex][1] = -( m_meanPermCoefficient * dHarmonic[1] + (1 - m_meanPermCoefficient) * dArithmetic[1] );

      connectionIndex++;
    }
  }


}


GEOSX_HOST_DEVICE
inline void SurfaceElementStencilWrapper::computeWeights( localIndex iconn,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar1,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar2,
                                                          real64 (& weight)[MAX_NUM_OF_CONNECTIONS][2],
                                                          real64 (& dWeight_dVar1 )[MAX_NUM_OF_CONNECTIONS][2],
                                                          real64 (& dWeight_dVar2 )[MAX_NUM_OF_CONNECTIONS][2] ) const
{
  // TODO: this should become star-delta method
  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k<numPointsInFlux( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    sumOfTrans += coefficient[er][esr][ei][0][0] * m_weights[iconn][k];
  }

  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {
      localIndex const er0  =  m_elementRegionIndices[iconn][k[0]];
      localIndex const esr0 =  m_elementSubRegionIndices[iconn][k[0]];
      localIndex const ei0  =  m_elementIndices[iconn][k[0]];

      localIndex const er1  =  m_elementRegionIndices[iconn][k[1]];
      localIndex const esr1 =  m_elementSubRegionIndices[iconn][k[1]];
      localIndex const ei1  =  m_elementIndices[iconn][k[1]];

      real64 const t0 = m_weights[iconn][0] * coefficient[er0][esr0][ei0][0][0]; // this is a bit insane to access perm
      real64 const t1 = m_weights[iconn][1] * coefficient[er1][esr1][ei1][0][0];

      real64 const harmonicWeight   = t0*t1 / sumOfTrans;
      real64 const arithmeticWeight = 0.25 * (t0+t1);

      real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

      weight[connectionIndex][0] = value;
      weight[connectionIndex][1] = -value;

      real64 const dt0_dvar1 = m_weights[iconn][0] * dCoeff_dVar1[er0][esr0][ei0][0][0];
      real64 const dt1_dvar1 = m_weights[iconn][1] * dCoeff_dVar1[er1][esr1][ei1][0][0];

      real64 dHarmonic_dvar1[2];
      dHarmonic_dvar1[0] = ( dt0_dvar1 * t1 * sumOfTrans - dt0_dvar1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
      dHarmonic_dvar1[1] = ( dt0_dvar1 * t1 * sumOfTrans - dt0_dvar1 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );

      real64 dArithmetic_dvar1[2];
      dArithmetic_dvar1[0] = 0.25 * dt0_dvar1;
      dArithmetic_dvar1[1] = 0.25 * dt1_dvar1;

      dWeight_dVar1[connectionIndex][0] =    m_meanPermCoefficient * dHarmonic_dvar1[0] + (1 - m_meanPermCoefficient) * dArithmetic_dvar1[0];
      dWeight_dVar1[connectionIndex][1] = -( m_meanPermCoefficient * dHarmonic_dvar1[1] + (1 - m_meanPermCoefficient) * dArithmetic_dvar1[1] );

      real64 const dt0_dvar2 = m_weights[iconn][0] * dCoeff_dVar2[er0][esr0][ei0][0][0];
      real64 const dt1_dvar2 = m_weights[iconn][1] * dCoeff_dVar2[er1][esr1][ei1][0][0];

      real64 dHarmonic_dvar2[2];
      dHarmonic_dvar2[0] = ( dt0_dvar2 * t1 * sumOfTrans - dt0_dvar2 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );
      dHarmonic_dvar2[1] = ( t0 * dt1_dvar2 * sumOfTrans - dt1_dvar2 * t0 * t1 ) / ( sumOfTrans * sumOfTrans );

      real64 dArithmetic_dvar2[2];
      dArithmetic_dvar2[0] = 0.25 * dt0_dvar2;
      dArithmetic_dvar2[1] = 0.25 * dt1_dvar2;

      dWeight_dVar2[connectionIndex][0] =   ( m_meanPermCoefficient * dHarmonic_dvar2[0] + (1 - m_meanPermCoefficient) * dArithmetic_dvar2[0] );
      dWeight_dVar2[connectionIndex][1] = -( m_meanPermCoefficient * dHarmonic_dvar2[1] + (1 - m_meanPermCoefficient) * dArithmetic_dvar2[1] );

      connectionIndex++;
    }
  }
}

GEOSX_HOST_DEVICE
inline void SurfaceElementStencilWrapper::computeWeights( localIndex iconn,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficientMultiplier,
                                                          R1Tensor const & gravityVector,
                                                          real64 (& weight)[MAX_NUM_OF_CONNECTIONS][2] ) const
{
  // TODO: this should become star-delta method
  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k<numPointsInFlux( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    real64 const mult = ( LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( m_cellCenterToEdgeCenters[iconn][k], gravityVector ) ) > MULTIPLIER_THRESHOLD )
      ? coefficientMultiplier[er][esr][ei][0][1] : coefficientMultiplier[er][esr][ei][0][0];

    sumOfTrans += mult * coefficient[er][esr][ei][0][0] * m_weights[iconn][k];
  }


  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {
      localIndex const er0  =  m_elementRegionIndices[iconn][k[0]];
      localIndex const esr0 =  m_elementSubRegionIndices[iconn][k[0]];
      localIndex const ei0  =  m_elementIndices[iconn][k[0]];

      localIndex const er1  =  m_elementRegionIndices[iconn][k[1]];
      localIndex const esr1 =  m_elementSubRegionIndices[iconn][k[1]];
      localIndex const ei1  =  m_elementIndices[iconn][k[1]];

      real64 const mult0 = ( LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( m_cellCenterToEdgeCenters[iconn][k[0]], gravityVector ) ) > MULTIPLIER_THRESHOLD )
  ? coefficientMultiplier[er0][esr0][ei0][0][1] : coefficientMultiplier[er0][esr0][ei0][0][0];
      real64 const mult1 = ( LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( m_cellCenterToEdgeCenters[iconn][k[1]], gravityVector ) ) > MULTIPLIER_THRESHOLD )
  ? coefficientMultiplier[er1][esr1][ei1][0][1] : coefficientMultiplier[er1][esr1][ei1][0][0];

      real64 const t0 = mult0 * m_weights[iconn][0] * coefficient[er0][esr0][ei0][0][0];
      real64 const t1 = mult1 * m_weights[iconn][1] * coefficient[er1][esr1][ei1][0][0];

      real64 const harmonicWeight   = t0*t1 / sumOfTrans;
      real64 const arithmeticWeight = 0.25 * (t0+t1);

      real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

      weight[connectionIndex][0] = value;
      weight[connectionIndex][1] = -value;

      connectionIndex++;
    }
  }
}

GEOSX_HOST_DEVICE
inline void SurfaceElementStencilWrapper::computeWeights( localIndex iconn,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficient1,
                                                          CoefficientAccessor< arrayView3d< real64 const > > const & coefficient1Multiplier,
                                                          CoefficientAccessor< arrayView1d< real64 const > > const & coefficient2,
                                                          R1Tensor const & unitGravityVector,
                                                          real64 ( & weight1 )[NUM_POINT_IN_FLUX],
                                                          real64 ( & weight2 )[NUM_POINT_IN_FLUX],
                                                          real64 ( & geometricWeight )[NUM_POINT_IN_FLUX] ) const
{
  real64 sumOfGeometricWeights = 0.0;

  for( localIndex k = 0; k < numPointsInFlux( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    real64 const cellToEdgeDistance = LvArray::tensorOps::l2Norm< 3 >( m_cellCenterToEdgeCenters[iconn][k] );
    real64 const edgeLength = m_weights[iconn][k] * cellToEdgeDistance;
    real64 const edgeToFaceDownDistance = -LvArray::tensorOps::AiBi< 3 >( m_cellCenterToEdgeCenters[iconn][k], unitGravityVector )
                                          * edgeLength / cellToEdgeDistance;

    real64 const mult = ( LvArray::math::abs( edgeToFaceDownDistance ) > MULTIPLIER_THRESHOLD )
      ? coefficient1Multiplier[er][esr][ei][0][1] : coefficient1Multiplier[er][esr][ei][0][0];

    weight1[k] = mult * coefficient1[er][esr][ei][0][0] * m_weights[iconn][k];
    weight2[k] = coefficient2[er][esr][ei] * edgeToFaceDownDistance;

    geometricWeight[k] = m_weights[iconn][k] / 12.0;
    sumOfGeometricWeights += geometricWeight[k];
  }

  for( localIndex k = 0; k < numPointsInFlux( iconn ); ++k )
  {
    geometricWeight[k] /= sumOfGeometricWeights;
  }
}


} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_ */
