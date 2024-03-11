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

#ifndef GEOS_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_
#define GEOS_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geos
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
 * @brief Describes properties of SurfaceElementStencil.
 *
 * This type of stencil allows for up to 6 surface elements to be connected in a flux computation.
 * The total number of pairwise connections is thus: 6*(6-1)/2 = 15.
 */
using SurfaceElementStencilTraits = StencilTraits< ArrayOfArrays, 6, 6, 15 >;

/**
 * @brief Provides access to the SurfaceElementStencil that may be called from a kernel function.
 */
class SurfaceElementStencilWrapper : public StencilWrapperBase< SurfaceElementStencilTraits >
{
public:

  static constexpr real64 avgWeights = 0.25;

  /// Threshold for the application of the permeability multiplier
  static constexpr real64 MULTIPLIER_THRESHOLD = 1e-10;

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
                                real64 const meanPermCoefficient );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex size() const { return m_elementRegionIndices.size(); }


  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOS_HOST_DEVICE
  localIndex stencilSize( localIndex index ) const override
  { return m_elementRegionIndices.sizeOfArray( index ); }



  using StencilWrapperBase< SurfaceElementStencilTraits >::computeWeights;
  /**
   * @brief Compute weights and derivatives w.r.t to one variable without coefficient
   * Used in ReactiveCompositionalMultiphaseOBL solver for thermal transmissibility computation:
   * here, conductivity is a part of operator and connot be used directly as a coefficient
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       real64 ( &weight )[maxNumConnections][2],
                       real64 ( &dWeight_dVar )[maxNumConnections][2] ) const;

  /**
   * @brief Compute weights and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] coefficientMultiplier view accessor to the coefficient multiplier used to compute the weights
   * @param[in] gravityVector gravity vector
   * @param[out] weight view weights
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const &coefficientMultiplier,
                       R1Tensor const &gravityVector,
                       real64 ( &weight )[maxNumConnections][2] ) const;

  /**
   * @brief Compute weights and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient1 view accessor to the first coefficient used to compute the first weights
   * @param[in] coefficient1Multiplier view accessor to the coefficient multiplier used to compute the first weights
   * @param[in] coefficient2 view accessor to the first coefficient used to compute the second weights
   * @param[in] gravityVector gravity vector
   * @param[out] weight1 view on the first weights
   * @param[out] weight2 view on the second weights
   * @param[out] geometricWeight view on the purely geometric weights
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const &coefficient1,
                       CoefficientAccessor< arrayView3d< real64 const > > const &coefficient1Multiplier,
                       CoefficientAccessor< arrayView1d< real64 const > > const &coefficient2,
                       R1Tensor const &gravityVector,
                       real64 ( &weight1 )[maxNumPointsInFlux],
                       real64 ( &weight2 )[maxNumPointsInFlux],
                       real64 ( &geometricWeight )[maxNumPointsInFlux] ) const;

  /**
   * @brief Compute the stabilization weights
   * @param[in] iconn connection index
   * @param[out] stabilizationWeight view weights
   */
  GEOS_HOST_DEVICE
  void computeStabilizationWeights( localIndex iconn,
                                    real64 ( & stabilizationWeight )[maxNumConnections][2] ) const
  {
    GEOS_UNUSED_VAR( iconn, stabilizationWeight );
  }

  /**
   * @brief Accessor to the CellCenterToEdgeCenter vector
   * @return the view const to the CellCenterToEdgeCenter vector
   */
  ArrayOfArraysView< R1Tensor const >
  getCellCenterToEdgeCenters() const { return m_cellCenterToEdgeCenters.toViewConst(); }

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

  /// Cell center to Edge center vector
  ArrayOfArraysView< R1Tensor > m_cellCenterToEdgeCenters;



  GEOS_HOST_DEVICE
  void
  computeWeightsBase( localIndex const iconn,
                      localIndex const (&k)[2],
                      localIndex const ielem,
                      arraySlice3d< real64 const > const & coefficient,
                      arraySlice3d< real64 const > const & dCoeff_dVar,
                      real64 & halfWeight,
                      real64 & dHalfWeight_dVar ) const override;

  GEOS_HOST_DEVICE
  void
  computeWeightsBase( localIndex const iconn,
                      localIndex const (&k)[2],
                      localIndex const ielem,
                      arraySlice3d< const real64 > const & coefficient,
                      arraySlice3d< const real64 > const & dCoeff_dVar1,
                      arraySlice3d< const real64 > const & dCoeff_dVar2,
                      real64 & halfWeight,
                      real64 ( & dHalfWeight_dVar )[2] ) const override
  {
    GEOS_UNUSED_VAR( iconn, k, ielem, coefficient, dCoeff_dVar1, dCoeff_dVar2, halfWeight, dHalfWeight_dVar );
  }

  GEOS_HOST_DEVICE
  void
    computeWeightsBase( localIndex const iconn,
                        localIndex const (&k)[2],
                        localIndex const ielem,
                        arraySlice3d< const real64 > const & coefficient,
                        arraySlice3d< const real64 > const & dCoeff_dVar1,
                        arraySlice4d< const real64 > const & dCoeff_dVar2,
                        real64 &halfWeight,
                        real64 ( &dHalfWeight_dVar )[2] ) const override;
};

/**
 * @brief Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class SurfaceElementStencil final : public StencilBase< SurfaceElementStencilTraits, SurfaceElementStencil >
{
public:

  virtual void move( LvArray::MemorySpace const space ) override;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

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
  using KernelWrapper = SurfaceElementStencilWrapper;

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
  real64 m_meanPermCoefficient = 1.0;

};

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::computeWeightsBase( const geos::localIndex iconn,
                                                  const geos::localIndex (& k)[2],
                                                  const localIndex ielem,
                                                  const arraySlice3d< const geos::real64 > & coefficient,
                                                  const arraySlice3d< const geos::real64 > & dCoeff_dVar,
                                                  geos::real64 & halfWeight,
                                                  geos::real64 & dHalfWeight_dVar ) const
{

  localIndex const ei  =  m_elementIndices[iconn][k[ielem]];
  halfWeight = m_weights[iconn][ielem] * coefficient[ei][0][0];
  dHalfWeight_dVar = m_weights[iconn][ielem] * dCoeff_dVar[ei][0][0];

}

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  computeWeights( localIndex iconn,
                  real64 ( & weight )[maxNumConnections][2],
                  real64 ( & dWeight_dVar )[maxNumConnections][2] ) const
{

  real64 sumOfTrans = 0.0;
  for( localIndex k=0; k<numPointsInFlux( iconn ); ++k )
  {
    sumOfTrans += m_weights[iconn][k];
  }

  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numPointsInFlux( iconn ); ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numPointsInFlux( iconn ); ++k[1] )
    {
      real64 const t0 = m_weights[iconn][0];
      real64 const t1 = m_weights[iconn][1];

      real64 const harmonicWeight   = t0*t1 / sumOfTrans;
      real64 const arithmeticWeight = 0.25 * (t0+t1);

      real64 const value = m_meanPermCoefficient * harmonicWeight + (1 - m_meanPermCoefficient) * arithmeticWeight;

      weight[connectionIndex][0] = value;
      weight[connectionIndex][1] = -value;

      real64 const dt0 = m_weights[iconn][0];
      real64 const dt1 = m_weights[iconn][1];

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



GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  computeWeightsBase( const geos::localIndex iconn,
                      const geos::localIndex (& k)[2],
                      const geos::localIndex ielem,
                      const arraySlice3d< const geos::real64 > & coefficient,
                      const arraySlice3d< const geos::real64 > & dCoeff_dVar1,
                      const arraySlice4d< const geos::real64 > & dCoeff_dVar2,
                      geos::real64 & halfWeight,
                      geos::real64 (& dHalfWeight_dVar) [2] ) const
{

  localIndex const ei  =  m_elementIndices[iconn][k[ielem]];
  halfWeight =
    m_weights[iconn][ielem] * coefficient[ei][0][0];         // this is a bit insane to access perm

  dHalfWeight_dVar[0] = m_weights[iconn][ielem] * dCoeff_dVar1[ei][0][0];
  dHalfWeight_dVar[1] = m_weights[iconn][ielem] * dCoeff_dVar2[ei][0][0][0];
}

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  computeWeights( localIndex iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficientMultiplier,
                  R1Tensor const & gravityVector,
                  real64 (& weight)[maxNumConnections][2] ) const
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



      connectionIndex++;
    }
  }
}

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  computeWeights( localIndex iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient1,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient1Multiplier,
                  CoefficientAccessor< arrayView1d< real64 const > > const & coefficient2,
                  R1Tensor const & unitGravityVector,
                  real64 ( & weight1 )[maxNumPointsInFlux],
                  real64 ( & weight2 )[maxNumPointsInFlux],
                  real64 ( & geometricWeight )[maxNumPointsInFlux] ) const
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

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  removeHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  for( localIndex k = 0; k < stencilSize( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    m_weights[iconn][k] = m_weights[iconn][k] / hydraulicAperture[er][esr][ei];
  }
}

GEOS_HOST_DEVICE
inline void
SurfaceElementStencilWrapper::
  addHydraulicApertureContribution( localIndex const iconn, ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > hydraulicAperture ) const
{
  for( localIndex k = 0; k < stencilSize( iconn ); ++k )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][k];
    localIndex const esr =  m_elementSubRegionIndices[iconn][k];
    localIndex const ei  =  m_elementIndices[iconn][k];

    m_weights[iconn][k] = m_weights[iconn][k] * hydraulicAperture[er][esr][ei];
  }
}

} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_ */
