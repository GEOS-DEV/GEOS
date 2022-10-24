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
 * @file SolidMechanicsThermoPoroElasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSTHERMOPOROELASTIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSTHERMOPOROELASTIC_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsExtrinsicData.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### ThermoPoroElastic Description
 * Implements the KernelBase interface functions required for using the
 * effective stress for the integration of the stress divergence. This is
 * templated on one of the "finite element kernel application" functions
 * such as geosx::finiteElement::RegionBasedKernelApplication.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermoPoroElastic :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  ThermoPoroElastic( NodeManager const & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     localIndex const targetRegionIndex,
                     SUBREGION_TYPE const & elementSubRegion,
                     FE_TYPE const & finiteElementSpace,
                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                     arrayView1d< globalIndex const > const inputDofNumber,
                     globalIndex const rankOffset,
                     CRSMatrixView< real64, globalIndex const > const inputMatrix,
                     arrayView1d< real64 > const inputRhs,
                     real64 const (&inputGravityVector)[3] ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.getExtrinsicData< extrinsicMeshData::solidMechanics::totalDisplacement >() ),
    m_uhat( nodeManager.getExtrinsicData< extrinsicMeshData::solidMechanics::incrementalDisplacement >() ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_density( inputConstitutiveType.getDensity() ),
    m_initialFluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::initialPressure >() ),
    m_fluidPressure_n( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >() ),
    m_fluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >() ),
    m_initialTemperature( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::initialTemperature >() ),
    m_temperature_n( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature_n >() ),
    m_temperature( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature >() )
  {}

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
                                       xLocal(),
                                       u_local(),
                                       uhat_local(),
                                       constitutiveStiffness()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the ThermoPoroElastic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
#endif
        stack.u_local[ a ][i] = m_disp[ localNodeIndex ][i];
        stack.uhat_local[ a ][i] = m_uhat[ localNodeIndex ][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 N[numNodesPerElem];
    real64 dNdX[ numNodesPerElem ][ 3 ];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    real64 strainIncrement[6]{};
    real64 totalStress[6]{};
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    // Evaluate total stress and its derivatives
    m_constitutiveUpdate.smallStrainUpdateThermalSinglePhase( k,
                                                              q,
                                                              m_initialFluidPressure[k],
                                                              m_fluidPressure_n[k],
                                                              m_fluidPressure[k],
                                                              m_initialTemperature[k],
                                                              m_temperature_n[k], 
                                                              m_temperature[k],
                                                              strainIncrement,
                                                              totalStress,
                                                              stiffness );

    for( localIndex i=0; i<6; ++i )
    {
      totalStress[i] *= -detJxW;
    }

    // Here we consider the bodyForce is purely from the solid
    real64 const bodyForce[3] = { m_gravityVector[0] * m_density( k, q )* detJxW,
                                  m_gravityVector[1] * m_density( k, q )* detJxW,
                                  m_gravityVector[2] * m_density( k, q )* detJxW };

    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     totalStress,
                                     N,
                                     bodyForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    // TODO: Does this work if BTDB is non-symmetric?
    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof =
          LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );
      }
    }


    return maxForce;
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];

  /// The rank global densities
  arrayView2d< real64 const > const m_density;

  /// The rank-global initial fluid pressure array
  arrayView1d< real64 const > const m_initialFluidPressure;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > const m_fluidPressure_n;
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global initial temperature array
  arrayView1d< real64 const > const m_initialTemperature;

  /// The rank-global temperature arrays.
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_temperature;

};

/// The factory used to construct a ThermoPoroElastic kernel.
using ThermoPoroElasticFactory = finiteElement::KernelFactory< ThermoPoroElastic,
                                                               arrayView1d< globalIndex const > const,
                                                               globalIndex,
                                                               CRSMatrixView< real64, globalIndex const > const,
                                                               arrayView1d< real64 > const,
                                                               real64 const (&)[3] >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSTHERMOPOROELASTIC_HPP_
