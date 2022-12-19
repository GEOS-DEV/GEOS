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
 * @file MultiphasePoromechanics_impl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_IMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_IMPL_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geosx
{

namespace poromechanicsKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
MultiphasePoromechanics( NodeManager const & nodeManager,
                         EdgeManager const & edgeManager,
                         FaceManager const & faceManager,
                         localIndex const targetRegionIndex,
                         SUBREGION_TYPE const & elementSubRegion,
                         FE_TYPE const & finiteElementSpace,
                         CONSTITUTIVE_TYPE & inputConstitutiveType,
                         arrayView1d< globalIndex const > const inputDispDofNumber,
                         globalIndex const rankOffset,
                         CRSMatrixView< real64, globalIndex const > const inputMatrix,
                         arrayView1d< real64 > const inputRhs,
                         real64 const (&inputGravityVector)[3],
                         string const inputFlowDofKey,
                         localIndex const numComponents,
                         localIndex const numPhases,
                         string const fluidModelKey ):
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
  m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_uhat( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
  m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
  m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
  m_solidDensity( inputConstitutiveType.getDensity() ),
  m_numComponents( numComponents ),
  m_numPhases( numPhases )
{
  GEOSX_ERROR_IF_GT_MSG( m_numComponents, numMaxComponents,
                         "MultiphasePoromechanics solver allows at most " << numMaxComponents << " components at the moment" );

  m_flowDofNumber = elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey );

  // extract fluid constitutive data views
  {
    string const fluidModelName = elementSubRegion.template getReference< string >( fluidModelKey );
    constitutive::MultiFluidBase const & fluid =
      elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelName );

    m_fluidPhaseDensity = fluid.phaseDensity();
    m_fluidPhaseDensity_n = fluid.phaseDensity_n();
    m_dFluidPhaseDensity = fluid.dPhaseDensity();

    m_fluidPhaseCompFrac = fluid.phaseCompFraction();
    m_fluidPhaseCompFrac_n = fluid.phaseCompFraction_n();
    m_dFluidPhaseCompFrac = fluid.dPhaseCompFraction();

    m_fluidPhaseMassDensity = fluid.phaseMassDensity();
    m_dFluidPhaseMassDensity = fluid.dPhaseMassDensity();

  }

  // extract views into flow solver data
  {
    using namespace fields::flow;

    m_fluidPressure_n = elementSubRegion.template getField< pressure_n >();
    m_fluidPressure = elementSubRegion.template getField< pressure >();

    m_fluidPhaseSaturation_n = elementSubRegion.template getField< phaseVolumeFraction_n >();

    m_fluidPhaseSaturation = elementSubRegion.template getField< phaseVolumeFraction >();
    m_dFluidPhaseSaturation = elementSubRegion.template getField< dPhaseVolumeFraction >();

    m_dGlobalCompFraction_dGlobalCompDensity =
      elementSubRegion.template getField< dGlobalCompFraction_dGlobalCompDensity >();
  }
}


/**
 * @brief Copy global values from primary field to a local stack array.
 * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
 *
 * For the MultiphasePoromechanics implementation, global values from the displacement,
 * incremental displacement, and degree of freedom numbers are placed into
 * element local stack storage.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  for( localIndex a=0; a<numSupportPoints; ++a )
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

  // Add stabilization to block diagonal parts of the local dResidualMomentum_dDisplacement (this
  // is a no-operation with FEM classes)
  real64 const stabilizationScaling = computeStabilizationScaling( k );
  m_finiteElementSpace.template addGradGradStabilizationMatrix
  < FE_TYPE, numDofPerTrialSupportPoint, false >( stack.feStack,
                                                  stack.dLocalResidualMomentum_dDisplacement,
                                                  -stabilizationScaling );
  m_finiteElementSpace.template
  addEvaluatedGradGradStabilizationVector< FE_TYPE,
                                           numDofPerTrialSupportPoint >
    ( stack.feStack,
    stack.uhat_local,
    reinterpret_cast< real64 (&)[numNodesPerElem][numDofPerTestSupportPoint] >(stack.localResidualMomentum),
    -stabilizationScaling );

}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{
  using namespace PDEUtilities;

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  real64 strainIncrement[6]{};
  real64 totalStress[6]{};
  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
  real64 dTotalStress_dPressure[6]{};
  real64 bodyForce[3]{};
  real64 dBodyForce_dVolStrainIncrement[3]{};
  real64 dBodyForce_dPressure[3]{};
  real64 dBodyForce_dComponents[3][numMaxComponents]{};
  real64 componentMassContentIncrement[numMaxComponents]{};
  real64 dComponentMassContent_dVolStrainIncrement[numMaxComponents]{};
  real64 dComponentMassContent_dPressure[numMaxComponents]{};
  real64 dComponentMassContent_dComponents[numMaxComponents][numMaxComponents]{};
  real64 poreVolumeConstraint;
  real64 dPoreVolumeConstraint_dPressure;
  real64 dPoreVolumeConstraint_dComponents[1][numMaxComponents]{};

  // Displacement finite element basis functions (N), basis function derivatives (dNdX), and
  // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
  real64 N[numNodesPerElem];
  real64 dNdX[numNodesPerElem][3];
  FE_TYPE::calcN( q, stack.feStack, N );
  real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                           stack.feStack, dNdX );

  // Compute strain increment
  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

  // Evaluate conserved quantities (total stress and fluid mass content) and their derivatives
  m_constitutiveUpdate.smallStrainUpdateMultiphase( k,
                                                    q,
                                                    NP,
                                                    NC,
                                                    m_fluidPressure_n[k],
                                                    m_fluidPressure[k],
                                                    strainIncrement,
                                                    m_gravityAcceleration,
                                                    m_gravityVector,
                                                    m_solidDensity( k, q ),
                                                    m_fluidPhaseDensity[k][q],
                                                    m_fluidPhaseDensity_n[k][q],
                                                    m_dFluidPhaseDensity[k][q],
                                                    m_fluidPhaseCompFrac[k][q],
                                                    m_fluidPhaseCompFrac_n[k][q],
                                                    m_dFluidPhaseCompFrac[k][q],
                                                    m_fluidPhaseMassDensity[k][q],
                                                    m_dFluidPhaseMassDensity[k][q],
                                                    m_fluidPhaseSaturation[k],
                                                    m_fluidPhaseSaturation_n[k],
                                                    m_dFluidPhaseSaturation[k],
                                                    m_dGlobalCompFraction_dGlobalCompDensity[k],
                                                    totalStress,
                                                    dTotalStress_dPressure,
                                                    bodyForce,
                                                    dBodyForce_dVolStrainIncrement,
                                                    dBodyForce_dPressure,
                                                    dBodyForce_dComponents,
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
    stiffness,   // fourth-order tensor handled via DiscretizationOps
    dNdX,
    -detJxW );

  if( m_gravityAcceleration > 0.0 )
  {
    BilinearFormUtilities::compute< displacementTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Divergence >
    (
      stack.dLocalResidualMomentum_dDisplacement,
      N,
      dBodyForce_dVolStrainIncrement,
      dNdX,
      detJxW );
  }

  // Compute local linear momentum balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< displacementTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::SymmetricGradient,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dPressure,
    dNdX,
    dTotalStress_dPressure,
    1.0,
    -detJxW );

  if( m_gravityAcceleration > 0.0 )
  {
    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dPressure,
      N,
      dBodyForce_dPressure,
      1.0,
      detJxW );
  }

  // Compute local linear momentum balance residual derivatives with respect to components
  if( m_gravityAcceleration > 0.0 )
  {
    BilinearFormUtilities::compute< displacementTestSpace,
                                    FunctionSpace::P0,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dComponents,
      N,
      dBodyForce_dComponents,
      1.0,
      detJxW );
  }

  // --- Mass balance equations
  // --- --- Local component mass balance residual
  LinearFormUtilities::compute< FunctionSpace::P0,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMass,
    1.0,
    componentMassContentIncrement,
    detJxW );

  // --- --- Compute local mass balance residual derivatives with respect to displacement
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  displacementTestSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Divergence >
  (
    stack.dLocalResidualMass_dDisplacement,
    1.0,
    dComponentMassContent_dVolStrainIncrement,
    dNdX,
    detJxW );

  // --- --- Compute local mass balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dPressure,
    1.0,
    dComponentMassContent_dPressure,
    1.0,
    detJxW );

  // --- --- Compute local mass balance residual derivatives with respect to components
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dComponents,
    1.0,
    dComponentMassContent_dComponents,
    1.0,
    detJxW );

  // --- Pore volume constraint equation
  // --- --- Local pore volume contraint residual
  LinearFormUtilities::compute< FunctionSpace::P0,
                                DifferentialOperator::Identity >
  (
    stack.localPoreVolumeConstraint,
    1.0,
    poreVolumeConstraint,
    detJxW );

  // --- --- Compute local pore volume contraint residual derivatives with respect to pressure
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalPoreVolumeConstraint_dPressure,
    1.0,
    dPoreVolumeConstraint_dPressure,
    1.0,
    detJxW );

  // --- --- Compute local pore volume contraint residual derivatives with respect to components
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalPoreVolumeConstraint_dComponents,
    1.0,
    dPoreVolumeConstraint_dComponents,
    1.0,
    detJxW );

}

/**
 * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  using namespace compositionalMultiphaseUtilities;

  GEOSX_UNUSED_VAR( k );

  real64 maxForce = 0;
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  int nUDof = numSupportPoints * numDofPerTestSupportPoint;
  constexpr int nMaxUDof = FE_TYPE::maxSupportPoints * numDofPerTestSupportPoint;

  // Apply equation/variable change transformation(s)
  real64 work[nMaxUDof > ( numMaxComponents + 1 ) ? nMaxUDof : numMaxComponents + 1];
  shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, nUDof, stack.dLocalResidualMass_dDisplacement, work );
  shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, 1, stack.dLocalResidualMass_dPressure, work );
  shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, m_numComponents, stack.dLocalResidualMass_dComponents, work );
  shiftElementsAheadByOneAndReplaceFirstElementWithSum( m_numComponents, stack.localResidualMass );

  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint * localNode + dim] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.dLocalResidualMomentum_dDisplacement[numDofPerTestSupportPoint * localNode + dim],
                                                                              nUDof );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] );
      maxForce = fmax( maxForce, fabs( stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] ) );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localPressureDofIndex,
                                                                              stack.dLocalResidualMomentum_dPressure[numDofPerTestSupportPoint * localNode + dim],
                                                                              1 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localComponentDofIndices,
                                                                              stack.dLocalResidualMomentum_dComponents[numDofPerTestSupportPoint * localNode + dim],
                                                                              m_numComponents );
    }
  }

  localIndex const dof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex[0] - m_dofRankOffset );
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

    m_matrix.template addToRow< serialAtomic >( dof + m_numComponents,
                                                stack.localPressureDofIndex,
                                                stack.dLocalPoreVolumeConstraint_dPressure[0],
                                                1 );

    m_matrix.template addToRow< serialAtomic >( dof + m_numComponents,
                                                stack.localComponentDofIndices,
                                                stack.dLocalPoreVolumeConstraint_dComponents[0],
                                                m_numComponents );

    RAJA::atomicAdd< serialAtomic >( &m_rhs[dof+m_numComponents], stack.localPoreVolumeConstraint[0] );
  }

  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
real64 MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeStabilizationScaling( localIndex const k ) const
{
  // TODO: generalize this to other constitutive models (currently we assume linear elasticity).
  return 2.0 * m_constitutiveUpdate.getShearModulus( k );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64 MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  GEOSX_MARK_FUNCTION;

  // Define a RAJA reduction variable to get the maximum residual contribution.
  RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

  forAll< POLICY >( numElems,
                    [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
    {
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    maxResidual.max( kernelComponent.complete( k, stack ) );
  } );
  return maxResidual.get();
}



} // namespace poromechanicsKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_IMPL_HPP_
