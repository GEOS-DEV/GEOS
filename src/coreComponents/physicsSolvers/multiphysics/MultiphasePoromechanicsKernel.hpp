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
 * @file MultiphasePoroelasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSKERNEL_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace PoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static multiphase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### MultiphasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static multiphase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class Multiphase :
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

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  static constexpr int numMaxComponents = 3;
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
  Multiphase( NodeManager const & nodeManager,
              EdgeManager const & edgeManager,
              FaceManager const & faceManager,
              localIndex const targetRegionIndex,
              SUBREGION_TYPE const & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType,
              arrayView1d< globalIndex const > const inputDispDofNumber,
              string const inputFlowDofKey,
              globalIndex const rankOffset,
              real64 const (&inputGravityVector)[3],
              localIndex const numComponents,
              localIndex const numPhases,
              string const fluidModelKey,
              CRSMatrixView< real64, globalIndex const > const inputMatrix,
              arrayView1d< real64 > const inputRhs ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition() ),
    m_disp( nodeManager.totalDisplacement() ),
    m_uhat( nodeManager.incrementalDisplacement() ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_numComponents( numComponents ),
    m_numPhases( numPhases )
  {
    GEOSX_ERROR_IF_GT_MSG( m_numComponents, numMaxComponents,
                           "MultiphasePoroelastic solver allows at most " << numMaxComponents << " components at the moment" );

    m_flowDofNumber = elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey );

    // extract fluid constitutive data views
    {
      string const fluidModelName = elementSubRegion.template getReference< string >( fluidModelKey );
      constitutive::MultiFluidBase const & fluid =
        elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelName );

      m_fluidPhaseDensity = fluid.phaseDensity();
      m_dFluidPhaseDensity_dPressure = fluid.dPhaseDensity_dPressure();
      m_dFluidPhaseDensity_dGlobalCompFraction = fluid.dPhaseDensity_dGlobalCompFraction();

      m_fluidPhaseCompFrac = fluid.phaseCompFraction();
      m_dFluidPhaseCompFrac_dPressure = fluid.dPhaseCompFraction_dPressure();
      m_dFluidPhaseCompFraction_dGlobalCompFraction = fluid.dPhaseCompFraction_dGlobalCompFraction();

      m_fluidPhaseMassDensity = fluid.phaseMassDensity();
      m_initialFluidTotalMassDensity = fluid.initialTotalMassDensity();

    }

    // extract views into flow solver data
    {
      using namespace extrinsicMeshData::flow;

      m_initialFluidPressure = elementSubRegion.template getExtrinsicData< initialPressure >();
      m_fluidPressureOld = elementSubRegion.template getExtrinsicData< pressure >();
      m_deltaFluidPressure = elementSubRegion.template getExtrinsicData< deltaPressure >();

      m_fluidPhaseDensityOld = elementSubRegion.template getExtrinsicData< phaseDensityOld >();
      m_fluidPhaseCompFracOld = elementSubRegion.template getExtrinsicData< phaseComponentFractionOld >();
      m_fluidPhaseSaturationOld = elementSubRegion.template getExtrinsicData< phaseVolumeFractionOld >();

      m_fluidPhaseSaturation = elementSubRegion.template getExtrinsicData< phaseVolumeFraction >();
      m_dFluidPhaseSaturation_dPressure = elementSubRegion.template getExtrinsicData< dPhaseVolumeFraction_dPressure >();
      m_dFluidPhaseSaturation_dGlobalCompDensity = elementSubRegion.template getExtrinsicData< dPhaseVolumeFraction_dGlobalCompDensity >();

      m_dGlobalCompFraction_dGlobalCompDensity =
        elementSubRegion.template getExtrinsicData< dGlobalCompFraction_dGlobalCompDensity >();
    }
  }

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

    static constexpr int numDispDofPerElem =  Base::StackVariables::numRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localResidualMomentum( Base::StackVariables::localResidual ),
      dLocalResidualMomentum_dDisplacement( Base::StackVariables::localJacobian ),
      dLocalResidualMomentum_dPressure{ {0.0} },
      localResidualMass{ 0.0 },
      dLocalResidualMass_dDisplacement{ {0.0} },
      dLocalResidualMass_dPressure{ {0.0} },
      dLocalResidualMass_dComponents{ {0.0} },
      localPoreVolumeConstraint{ 0.0 },
      dLocalPoreVolumeConstraint_dPressure{ {0.0} },
      dLocalPoreVolumeConstraint_dComponents{ {0.0} },
      localPressureDofIndex{ 0 },
      localComponentDofIndices{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 ( &localResidualMomentum )[numDispDofPerElem];
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1];

    real64 localResidualMass[numMaxComponents];
    real64 dLocalResidualMass_dDisplacement[numMaxComponents][numDispDofPerElem];
    real64 dLocalResidualMass_dPressure[numMaxComponents][1];
    real64 dLocalResidualMass_dComponents[numMaxComponents][numMaxComponents];

    real64 localPoreVolumeConstraint[1];
    real64 dLocalPoreVolumeConstraint_dPressure[1][1];
    real64 dLocalPoreVolumeConstraint_dComponents[1][numMaxComponents];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localPressureDofIndex[1];
    globalIndex localComponentDofIndices[numMaxComponents];

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the Multiphase implementation, global values from the displacement,
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
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    stack.localPressureDofIndex[0] = m_flowDofNumber[k];
    for( int flowDofIndex=0; flowDofIndex < numMaxComponents; ++flowDofIndex )
    {
      stack.localComponentDofIndices[flowDofIndex] = stack.localPressureDofIndex[0] + flowDofIndex + 1;
    }

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    using namespace PDEUtilities;

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
//    constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;

    real64 strainIncrement[6]{};
    real64 totalStress[6]{};
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    real64 dTotalStress_dPressure[6]{};
    real64 bodyForce[3]{};
    real64 dBodyForce_dVolStrainIncrement[3]{};
    real64 dBodyForce_dPressure[3]{};
    real64 componentMassContentIncrement[numMaxComponents]{};
    real64 dComponentMassContent_dVolStrainIncrement[numMaxComponents]{};
    real64 dComponentMassContent_dPressure[numMaxComponents]{};
    real64 dComponentMassContent_dComponents[numMaxComponents][numMaxComponents]{};
    real64 poreVolumeConstraint;
    real64 dPoreVolumeConstraint_dPressure;
    real64 dPoreVolumeConstraint_dComponents[numMaxComponents]{};

    // Displacement finite element basis functions (N), basis function derivatives (dNdX), and
    // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Compute strain increment
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    // Evaluate conserved quantities (total stress and fluid mass content) and their derivatives
    m_constitutiveUpdate.smallStrainUpdateMultiphase( k,
                                                      q,
                                                      NP,
                                                      NC,
                                                      m_initialFluidPressure[k],
                                                      m_fluidPressureOld[k],
                                                      m_deltaFluidPressure[k],
                                                      strainIncrement,
                                                      m_gravityAcceleration,
                                                      m_gravityVector,
                                                      m_solidDensity( k, q ),
                                                      m_initialFluidTotalMassDensity( k, q ),
                                                      m_fluidPhaseDensity[k][q],
                                                      m_fluidPhaseDensityOld[k],
                                                      m_dFluidPhaseDensity_dPressure[k][q],
                                                      m_dFluidPhaseDensity_dGlobalCompFraction[k][q],
                                                      m_fluidPhaseCompFrac[k][q],
                                                      m_fluidPhaseCompFracOld[k],
                                                      m_dFluidPhaseCompFrac_dPressure[k][q],
                                                      m_dFluidPhaseCompFraction_dGlobalCompFraction[k][q],
                                                      m_fluidPhaseMassDensity[k][q],
                                                      m_fluidPhaseSaturation[k],
                                                      m_fluidPhaseSaturationOld[k],
                                                      m_dFluidPhaseSaturation_dPressure[k],
                                                      m_dFluidPhaseSaturation_dGlobalCompDensity[k],
                                                      m_dGlobalCompFraction_dGlobalCompDensity[k],
                                                      totalStress,
                                                      dTotalStress_dPressure,
                                                      bodyForce,
                                                      dBodyForce_dVolStrainIncrement,
                                                      dBodyForce_dPressure,
                                                      componentMassContentIncrement,
                                                      dComponentMassContent_dVolStrainIncrement,
                                                      dComponentMassContent_dPressure,
                                                      dComponentMassContent_dComponents,
                                                      stiffness,
                                                      poreVolumeConstraint,
                                                      dPoreVolumeConstraint_dPressure,
                                                      dPoreVolumeConstraint_dComponents );

    // Compute local linear momentum balance residual
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::SymmetricGradient >
    (
      stack.localResidualMomentum,
      dNdX,
      totalStress,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      LinearFormUtilities::compute< displacementTestSpace,
                                    DifferentialOperator::Identity >
      (
        stack.localResidualMomentum,
        N,
        bodyForce,
        detJxW );
    }

    // Compute local linear momentum balance residual derivatives with respect to displacement
    BilinearFormUtilities::compute< displacementTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::SymmetricGradient >
    (
      stack.dLocalResidualMomentum_dDisplacement,
      dNdX,
      stiffness, // fourth-order tensor handled via DiscretizationOps
      dNdX,
      -detJxW );

    // --- dBodyForce_dVoldStrain derivative neglected

    // Compute local linear momentum balance residual derivatives with respect to pressure
    real64 Np[1] = { 1.0 };

    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dPressure,
      dNdX,
      dTotalStress_dPressure,
      Np,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      // dBodyForce_dVar neglected. with Var = {dVolStrain, dPressure, dComponents}
      // Assumptions: (  i) dMixtureDens_dVolStrain contribution is neglected
      //              ( ii) grains are assumed incompressible
      //              (iii) TODO add dMixtureDens_dPressure and dMixtureDens_dGlobalCompDensity
    }

    // --- Mass balance equations
    // --- --- loop over components
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      // Local component mass balance residual
      stack.localResidualMass[ic] = stack.localResidualMass[ic]
                                    + componentMassContentIncrement[ic] * detJxW;

      // Compute local mass balance residual derivatives with respect to displacement
      for( integer a = 0; a < numNodesPerElem; ++a )
      {
        stack.dLocalResidualMass_dDisplacement[ic][a*3+0] = stack.dLocalResidualMass_dDisplacement[ic][a*3+0]
                                                            + dNdX[a][0] * dComponentMassContent_dVolStrainIncrement[ic] * detJxW;
        stack.dLocalResidualMass_dDisplacement[ic][a*3+1] = stack.dLocalResidualMass_dDisplacement[ic][a*3+1]
                                                            + dNdX[a][1] * dComponentMassContent_dVolStrainIncrement[ic] * detJxW;
        stack.dLocalResidualMass_dDisplacement[ic][a*3+2] = stack.dLocalResidualMass_dDisplacement[ic][a*3+2]
                                                            + dNdX[a][2] * dComponentMassContent_dVolStrainIncrement[ic] * detJxW;
      }

      // Compute local mass balance residual derivatives with respect to pressure
      stack.dLocalResidualMass_dPressure[ic][0] = stack.dLocalResidualMass_dPressure[ic][0]
                                                  + dComponentMassContent_dPressure[ic] * detJxW;

      // Compute local mass balance residual derivatives with respect to components
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        stack.dLocalResidualMass_dComponents[ic][jc] = stack.dLocalResidualMass_dComponents[ic][jc]
                                                       + dComponentMassContent_dComponents[ic][jc] * detJxW;
      }
    }

    // --- Volume balance equation
    // sum contributions to component accumulation from each phase
    stack.localPoreVolumeConstraint[0] = stack.localPoreVolumeConstraint[0]
                                         + poreVolumeConstraint * detJxW;
    stack.dLocalPoreVolumeConstraint_dPressure[0][0] = stack.dLocalPoreVolumeConstraint_dPressure[0][0]
                                                       + dPoreVolumeConstraint_dPressure * detJxW;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      stack.dLocalPoreVolumeConstraint_dComponents[0][jc] = stack.dLocalPoreVolumeConstraint_dComponents[0][jc]
                                                            + dPoreVolumeConstraint_dComponents[jc] * detJxW;
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    using namespace CompositionalMultiphaseUtilities;

    GEOSX_UNUSED_VAR( k );

    real64 maxForce = 0;

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.dLocalResidualMomentum_dDisplacement );

    //int nFlowDof = m_numComponents + 1;
    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    // Apply equation/variable change transformation(s)
    real64 work[nUDof > ( numMaxComponents + 1 ) ? nUDof : numMaxComponents + 1];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, nUDof, stack.dLocalResidualMass_dDisplacement, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, 1, stack.dLocalResidualMass_dPressure, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, m_numComponents, stack.dLocalResidualMass_dComponents, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( m_numComponents, stack.localResidualMass );

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint * localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.dLocalResidualMomentum_dDisplacement[numDofPerTestSupportPoint * localNode + dim],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax( maxForce, fabs( stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] ) );

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localPressureDofIndex,
                                                                                stack.dLocalResidualMomentum_dPressure[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );

      }
    }

    localIndex dof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex[0] - m_dofRankOffset );
    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      for( localIndex i = 0; i < m_numComponents; ++i )
      {
        m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof + i,
                                                                        stack.localRowDofIndex,
                                                                        stack.dLocalResidualMass_dDisplacement[i],
                                                                        nUDof );
        m_matrix.template addToRow< serialAtomic >( dof + i,
                                                    stack.localPressureDofIndex,
                                                    stack.dLocalResidualMass_dPressure[i],
                                                    1 );
        m_matrix.template addToRow< serialAtomic >( dof + i,
                                                    stack.localComponentDofIndices,
                                                    stack.dLocalResidualMass_dComponents[i],
                                                    m_numComponents );
        RAJA::atomicAdd< serialAtomic >( &m_rhs[dof+i], stack.localResidualMass[i] );
      }
    }

    dof = dof + m_numComponents;
    if( 0 <= dof && dof < m_matrix.numRows() )
    {

      m_matrix.template addToRow< serialAtomic >( dof,
                                                  stack.localPressureDofIndex,
                                                  stack.dLocalPoreVolumeConstraint_dPressure[0],
                                                  1 );

      m_matrix.template addToRow< serialAtomic >( dof,
                                                  stack.localComponentDofIndices,
                                                  stack.dLocalPoreVolumeConstraint_dComponents[0],
                                                  m_numComponents );

      RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localPoreVolumeConstraint[0] );
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global density
  arrayView2d< real64 const > m_solidDensity;

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity;
  arrayView2d< real64 const, compflow::USD_PHASE > m_fluidPhaseDensityOld;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_dFluidPhaseDensity_dPressure;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseDensity_dGlobalCompFraction;

  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_COMP > m_fluidPhaseCompFracOld;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_dFluidPhaseCompFrac_dPressure;

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseMassDensity;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_dFluidPhaseMassDensity_dPressure;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseMassDensity_dGlobalCompFraction;

  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > m_initialFluidTotalMassDensity;

  arrayView2d< real64 const, compflow::USD_PHASE > m_fluidPhaseSaturation;
  arrayView2d< real64 const, compflow::USD_PHASE > m_fluidPhaseSaturationOld;
  arrayView2d< real64 const, compflow::USD_PHASE > m_dFluidPhaseSaturation_dPressure;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dFluidPhaseSaturation_dGlobalCompFraction;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dFluidPhaseSaturation_dGlobalCompDensity;

  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dGlobalCompFraction_dGlobalCompDensity;

  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > m_dFluidPhaseCompFraction_dGlobalCompFraction;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > m_flowDofNumber;

  /// The rank-global initial fluid pressure array.
  arrayView1d< real64 const > m_initialFluidPressure;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > m_fluidPressureOld;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > m_deltaFluidPressure;

  /// Number of components
  localIndex const m_numComponents;

  /// Number of phases
  localIndex const m_numPhases;

};

using MultiphaseKernelFactory = finiteElement::KernelFactory< Multiphase,
                                                              arrayView1d< globalIndex const > const,
                                                              string const,
                                                              globalIndex const,
                                                              real64 const (&)[3],
                                                              localIndex const,
                                                              localIndex const,
                                                              string const,
                                                              CRSMatrixView< real64, globalIndex const > const,
                                                              arrayView1d< real64 > const >;

} // namespace PoroelasticKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSKERNEL_HPP_
