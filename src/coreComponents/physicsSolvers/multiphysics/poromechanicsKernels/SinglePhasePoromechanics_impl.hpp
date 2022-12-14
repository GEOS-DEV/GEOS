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
 * @file SinglePhasePoromechanics_impl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_

#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geosx
{

namespace poromechanicsKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
SinglePhasePoromechanics( NodeManager const & nodeManager,
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
  m_X( nodeManager.referencePosition()),
  m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_uhat( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
  m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
  m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
  m_solidDensity( inputConstitutiveType.getDensity() ),
  m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
  m_fluidDensity_n( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
  m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).dDensity_dPressure() ),
  m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
  m_fluidPressure_n( elementSubRegion.template getField< fields::flow::pressure_n >() ),
  m_fluidPressure( elementSubRegion.template getField< fields::flow::pressure >() )
{}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
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

  stack.localFlowDofIndex[0] = m_flowDofNumber[k];

}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{
  // Governing equations (strong form)
  // ---------------------------------
  //
  //   divergence( totalStress ) + bodyForce = 0                       (quasi-static linear momentum balance)
  //   dFluidMassContent_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
  //
  // with currently the following dependencies on the strainIncrement tensor and pressure
  //
  //   totalStress      = totalStress( strainIncrement, pressure)
  //   bodyForce        = bodyForce( strainIncrement, pressure)
  //   fluidMassContent = fluidMassContent( strainIncrement, pressure)
  //   fluidMassFlux    = fludiMassFlux( pressure)
  //
  // Note that the fluidMassFlux will depend on the straiIncrement if a stress-dependent constitutive
  // model is assumed. A dependency on pressure can also occur in the source term of the mass
  // balance equation, e.g. if a Peaceman well model is used.
  //
  // In this kernel cell-based contributions to Jacobian matrix and residual
  // vector are computed. The face-based contributions, namely associated with the fluid
  // mass flux term are computed in a different kernel. The source term in the mass balance
  // equation is also treated elsewhere.
  //
  // Integration in time is performed using a backward Euler scheme (in the mass balance
  // equation LHS and RHS are multiplied by the timestep).
  //
  // Using a weak formulation of the governing equation the following terms are assembled in this kernel
  //
  //   Rmom = - \int symmetricGradient( \eta ) : totalStress + \int \eta \cdot bodyForce = 0
  //   Rmas = \int \chi ( fluidMassContent - fluidMassContent_n) = 0
  //
  //   dRmom_dVolStrain = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dVolStrain
  //                      + \int \eta \cdot dBodyForce_dVolStrain
  //   dRmom_dPressure  = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dPressure
  //                      + \int \eta \cdot dBodyForce_dPressure
  //   dRmas_dVolStrain = \int \chi dFluidMassContent_dVolStrain
  //   dRmas_dPressure  = \int \chi dFluidMassContent_dPressure
  //
  // with \eta and \chi test basis functions for the displacement and pressure field, respectively.
  // A continuous interpolation is used for the displacement, with \eta continuous finite element
  // basis functions. A piecewise-constant approximation is used for the pressure.

  using namespace PDEUtilities;

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
  constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

  real64 strainIncrement[6]{};
  real64 totalStress[6]{};
  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;   // Could this be called dTotalStress_dStrainIncrement?
  real64 dTotalStress_dPressure[6]{};
  real64 bodyForce[3]{};
  real64 dBodyForce_dVolStrainIncrement[3]{};
  real64 dBodyForce_dPressure[3]{};
  real64 fluidMassContentIncrement;
  real64 dFluidMassContent_dPressure;
  real64 dFluidMassContent_dVolStrainIncrement;

  // Displacement finite element basis functions (N), basis function derivatives (dNdX), and
  // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
  real64 N[numNodesPerElem];
  real64 dNdX[numNodesPerElem][3];
  FE_TYPE::calcN( q, N );
  real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

  // Compute strain increment
  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

  // Evaluate conserved quantities (total stress and fluid mass content) and their derivatives
  m_constitutiveUpdate.smallStrainUpdateSinglePhase( k,
                                                     q,
                                                     m_fluidPressure_n[k],
                                                     m_fluidPressure[k],
                                                     strainIncrement,
                                                     m_gravityAcceleration,
                                                     m_gravityVector,
                                                     m_solidDensity( k, q ),
                                                     m_fluidDensity_n( k, q ),
                                                     m_fluidDensity( k, q ),
                                                     m_dFluidDensity_dPressure( k, q ),
                                                     totalStress,
                                                     dTotalStress_dPressure,
                                                     bodyForce,
                                                     dBodyForce_dVolStrainIncrement,
                                                     dBodyForce_dPressure,
                                                     fluidMassContentIncrement,
                                                     dFluidMassContent_dPressure,
                                                     dFluidMassContent_dVolStrainIncrement,
                                                     stiffness );

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


  // Compute local mass balance residual
  LinearFormUtilities::compute< pressureTestSpace,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMass,
    1.0,
    fluidMassContentIncrement,
    detJxW );

  // Compute local mass balance residual derivatives with respect to displacement
  BilinearFormUtilities::compute< pressureTestSpace,
                                  displacementTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Divergence >
  (
    stack.dLocalResidualMass_dDisplacement,
    1.0,
    dFluidMassContent_dVolStrainIncrement,
    dNdX,
    detJxW );

  // Compute local mass balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dPressure,
    1.0,
    dFluidMassContent_dPressure,
    1.0,
    detJxW );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
//    constexpr FunctionSpace pressureTestSpace = pressureTrialSpace
  GEOSX_UNUSED_VAR( k );
  real64 maxForce = 0;

  constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

  for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.dLocalResidualMomentum_dDisplacement[numDofPerTestSupportPoint * localNode + dim],
                                                                              numNodesPerElem * numDofPerTrialSupportPoint );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] );
      maxForce = fmax( maxForce, fabs( stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] ) );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localFlowDofIndex,
                                                                              stack.dLocalResidualMomentum_dPressure[numDofPerTestSupportPoint * localNode + dim],
                                                                              1 );

    }
  }


  localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[0] - m_dofRankOffset );
  if( 0 <= dof && dof < m_matrix.numRows() )
  {
    m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof,
                                                                    stack.localRowDofIndex,
                                                                    stack.dLocalResidualMass_dDisplacement[0],
                                                                    nUDof );
    m_matrix.template addToRow< serialAtomic >( dof,
                                                stack.localFlowDofIndex,
                                                stack.dLocalResidualMass_dPressure[0],
                                                1 );
    RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localResidualMass[0] );
  }

  return maxForce;
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_
