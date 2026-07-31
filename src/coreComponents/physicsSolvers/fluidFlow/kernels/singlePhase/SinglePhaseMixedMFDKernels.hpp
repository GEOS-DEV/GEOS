/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseMixedMFDKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "finiteVolume/mimeticInnerProducts/AdaptiveInnerProduct.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "mixedMimetic/MixedMimeticDispatch.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/adaptivity/GlobalAdaptationKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{
namespace singlePhaseMixedMFDKernels
{

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_FACE number of faces per element
 * @tparam IP the type of MFD inner product used in the MFD-compatible cells
 * @brief Assemble the cell-wise contributions of the mixed mimetic saddle-point system:
 *        the face-based constitutive (Darcy closure) rows and the cell-based divergence terms.
 *
 * On each cell, the local mass matrix is selected according to the adaptive stencil flag:
 * TPFA-compatible cells use the diagonal TPFA mass matrix, MFD-compatible cells use IP.
 * The face rows are assembled with respect to the fixed global orientation of the face normal,
 * so that the (unknown) interface pressure cancels across interior interfaces.
 */
template< integer NUM_FACE, typename IP >
class ElementBasedAssemblyKernel
{
public:

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;

  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              real64 const & lengthTolerance,
                              string const elemDofKey,
                              string const faceDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              CellElementSubRegion const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::PermeabilityBase const & permeability,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_lengthTolerance( lengthTolerance ),
    m_dt( dt ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_elemDofNumber( subRegion.getReference< array1d< globalIndex > >( elemDofKey ) ),
    m_faceGhostRank( faceManager.ghostRank() ),
    m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
    m_elemToFaces( subRegion.faceList().toViewConst() ),
    m_elemCenter( subRegion.getElementCenter() ),
    m_elemVolume( subRegion.getElementVolume() ),
    m_elemGravCoef( subRegion.getField< fields::flow::gravityCoefficient >() ),
    m_faceToNodes( faceManager.nodeList().toViewConst() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_faceCenter( faceManager.faceCenter() ),
    m_faceNormal( faceManager.faceNormal() ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_nodePosition( nodeManager.referencePosition() ),
    m_elemPerm( permeability.permeability() ),
    m_elemPres( subRegion.getField< fields::flow::pressure >() ),
    m_faceFlux( faceManager.getField< fields::mixedMimetic::faceMassFlux >() ),
    m_bcPres( faceManager.getField< fields::flow::bcPressure >() ),
    m_isPresBcFace( faceManager.getField< fields::flow::isBoundaryFace >() ),
    m_stencilFlag( subRegion.getField< fields::mixedMimetic::stencilFlag >() ),
    m_elemDens( fluid.density() ),
    m_dElemDens( fluid.dDensity() ),
    m_mob( subRegion.getField< fields::flow::mobility >() ),
    m_dMob( subRegion.getField< fields::flow::dMobility >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
      : massMatrix( NUM_FACE, NUM_FACE )
    {}

    /// adaptive local mass matrix: chi * M_mfd + (1 - chi) * M_tpfa
    stackArray2d< real64, NUM_FACE *NUM_FACE > massMatrix;

    /// sign relating the cell-outward normal to the fixed global face normal
    real64 orientation[NUM_FACE]{};

    /// 1 if the face is a domain-boundary face without a pressure boundary condition (no-flow)
    integer isNoFlowFace[NUM_FACE]{};

    /// face constitutive residuals (local, outward-oriented) and their derivatives
    real64 faceResidual[NUM_FACE]{};
    real64 dFaceResidual_dFlux[NUM_FACE][NUM_FACE]{};
    real64 dFaceResidual_dPres[NUM_FACE]{};

    /// divergence of the face fluxes and its derivatives
    real64 divFlux = 0.0;
    real64 dDivFlux_dFlux[NUM_FACE]{};

    localIndex cellCenteredEqnRowIndex = 0;
    localIndex faceCenteredEqnRowIndex[NUM_FACE]{};
    globalIndex elemDofColIndex = 0;
    globalIndex faceDofColIndices[NUM_FACE]{};
  };

  /**
   * @brief Performs the setup phase for the kernel.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.cellCenteredEqnRowIndex = m_elemDofNumber[ei] - m_rankOffset;
    stack.elemDofColIndex = m_elemDofNumber[ei];
    for( integer iFaceLoc = 0; iFaceLoc < NUM_FACE; ++iFaceLoc )
    {
      localIndex const kf = m_elemToFaces[ei][iFaceLoc];
      stack.faceCenteredEqnRowIndex[iFaceLoc] = m_faceDofNumber[kf] - m_rankOffset;
      stack.faceDofColIndices[iFaceLoc] = m_faceDofNumber[kf];
    }
  }

  /**
   * @brief Compute the local constitutive rows and the flux divergence in a given element.
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    real64 const perm[ 3 ] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };

    // step 1: adaptive inner product, M = chi * M_mfd + (1 - chi) * M_tpfa
    real64 const chi = static_cast< real64 >( m_stencilFlag[ei] );
    mimeticInnerProduct::AdaptiveInnerProduct< IP >::template computeM< NUM_FACE >( m_nodePosition,
                                                                                    m_faceToNodes,
                                                                                    m_elemToFaces[ei],
                                                                                    m_elemCenter[ei],
                                                                                    m_elemVolume[ei],
                                                                                    perm,
                                                                                    m_lengthTolerance,
                                                                                    chi,
                                                                                    stack.massMatrix );

    // step 2: cell mobility scaling of the inverse Darcy coefficient
    real64 const mob = m_mob[ei];
    real64 const dMob_dPres = m_dMob[ei][DerivOffset::dP];
    real64 const invMob = 1.0 / mob;
    real64 const dInvMob_dPres = -dMob_dPres * invMob * invMob;

    real64 const ccDens = m_elemDens[ei][0];
    real64 const dCcDens_dPres = m_dElemDens[ei][0][DerivOffset::dP];

    // step 3: face orientations, localized fluxes and no-flow flags
    real64 localFlux[NUM_FACE]{};
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];

      real64 dotProd = 0.0;
      for( integer d = 0; d < 3; ++d )
      {
        dotProd += ( m_faceCenter[kf][d] - m_elemCenter[ei][d] ) * m_faceNormal[kf][d];
      }
      stack.orientation[i] = ( dotProd >= 0.0 ) ? 1.0 : -1.0;
      localFlux[i] = stack.orientation[i] * m_faceFlux[kf];

      bool const onBoundary = ( m_elemRegionList[kf][0] == -1 || m_elemRegionList[kf][1] == -1 );
      stack.isNoFlowFace[i] = ( onBoundary && m_isPresBcFace[kf] == 0 ) ? 1 : 0;
    }

    // step 4: constitutive (Darcy closure) rows:
    //   invMob * sum_j M_ij (sigma_j m_j) - p_K + pi_f + rho * ( g_K - g_f ) = 0
    // the interface pressure pi_f is dropped on interior faces (it cancels upon global
    // orientation-signed assembly) and replaced by the prescribed value on boundary faces
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      if( stack.isNoFlowFace[i] == 1 )
      {
        continue;
      }

      localIndex const kf = m_elemToFaces[ei][i];

      real64 mDotFlux = 0.0;
      for( integer j = 0; j < NUM_FACE; ++j )
      {
        mDotFlux += stack.massMatrix( i, j ) * localFlux[j];
        stack.dFaceResidual_dFlux[i][j] = invMob * stack.massMatrix( i, j ) * stack.orientation[j];
      }

      real64 const gravCoefDif = m_elemGravCoef[ei] - m_faceGravCoef[kf];

      stack.faceResidual[i] = invMob * mDotFlux - m_elemPres[ei] + ccDens * gravCoefDif;
      stack.dFaceResidual_dPres[i] = dInvMob_dPres * mDotFlux - 1.0 + dCcDens_dPres * gravCoefDif;

      bool const onBoundary = ( m_elemRegionList[kf][0] == -1 || m_elemRegionList[kf][1] == -1 );
      if( onBoundary && m_isPresBcFace[kf] == 1 )
      {
        stack.faceResidual[i] += m_bcPres[kf];
      }
    }

    // step 5: divergence of the face fluxes for the cell-centered mass conservation equation
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      stack.divFlux += m_dt * localFlux[i];
      stack.dDivFlux_dFlux[i] = m_dt * stack.orientation[i];
    }
  }

  /**
   * @brief Assemble the local contributions into the global system.
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // step 1: cell-centered mass conservation contribution (single thread per row)
    if( m_elemGhostRank[ei] < 0 )
    {
      m_localRhs[stack.cellCenteredEqnRowIndex] += stack.divFlux;

      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( stack.cellCenteredEqnRowIndex,
                                                                  &stack.faceDofColIndices[0],
                                                                  &stack.dDivFlux_dFlux[0],
                                                                  NUM_FACE );
    }

    // step 2: face-centered constitutive rows, signed by the global face orientation
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      if( m_faceGhostRank[kf] >= 0 )
      {
        continue;
      }

      if( stack.isNoFlowFace[i] == 1 )
      {
        // essential no-flow condition: m_f = 0 (row owned by the unique adjacent cell)
        real64 const one = 1.0;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[stack.faceCenteredEqnRowIndex[i]], m_faceFlux[kf] );
        m_localMatrix.addToRow< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[i],
                                                        &stack.faceDofColIndices[i],
                                                        &one,
                                                        1 );
        continue;
      }

      real64 const sigma = stack.orientation[i];

      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[stack.faceCenteredEqnRowIndex[i]], sigma * stack.faceResidual[i] );

      real64 const dRes_dPres = sigma * stack.dFaceResidual_dPres[i];
      m_localMatrix.addToRow< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[i],
                                                      &stack.elemDofColIndex,
                                                      &dRes_dPres,
                                                      1 );

      real64 dRes_dFlux[NUM_FACE]{};
      for( integer j = 0; j < NUM_FACE; ++j )
      {
        dRes_dFlux[j] = sigma * stack.dFaceResidual_dFlux[i][j];
      }
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[i],
                                                                          &stack.faceDofColIndices[0],
                                                                          &dRes_dFlux[0],
                                                                          NUM_FACE );
    }
  }

  /**
   * @brief Performs the kernel launch
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// offset for my MPI rank
  globalIndex const m_rankOffset;

  /// length tolerance
  real64 const m_lengthTolerance;

  /// time step size
  real64 const m_dt;

  /// ghost ranks and dof numbers
  arrayView1d< integer const > const m_elemGhostRank;
  arrayView1d< globalIndex const > const m_elemDofNumber;
  arrayView1d< integer const > const m_faceGhostRank;
  arrayView1d< globalIndex const > const m_faceDofNumber;

  /// topological and geometrical data
  arrayView2d< localIndex const > const m_elemToFaces;
  arrayView2d< real64 const > const m_elemCenter;
  arrayView1d< real64 const > const m_elemVolume;
  arrayView1d< real64 const > const m_elemGravCoef;
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  arrayView1d< real64 const > const m_faceGravCoef;
  arrayView2d< real64 const > const m_faceCenter;
  arrayView2d< real64 const > const m_faceNormal;
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;

  /// permeability
  arrayView3d< real64 const > const m_elemPerm;

  /// primary variables and boundary data
  arrayView1d< real64 const > const m_elemPres;
  arrayView1d< real64 const > const m_faceFlux;
  arrayView1d< real64 const > const m_bcPres;
  arrayView1d< integer const > const m_isPresBcFace;

  /// adaptive stencil flag
  arrayView1d< integer const > const m_stencilFlag;

  /// fluid data
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_elemDens;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dElemDens;
  arrayView1d< real64 const > const m_mob;
  arrayView2d< real64 const > const m_dMob;

  /// View on the local CRS matrix and RHS
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  arrayView1d< real64 > const m_localRhs;

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   CellElementSubRegion const & subRegion,
                   mimeticInnerProduct::MimeticInnerProductBase const & mimeticInnerProductBase,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::PermeabilityBase const & permeability,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    mixedMimeticInnerProductDispatch( mimeticInnerProductBase,
                                      [&] ( auto const mimeticInnerProduct )
    {
      using IP = TYPEOFREF( mimeticInnerProduct );

      mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        ElementBasedAssemblyKernel< NUM_FACES, IP >
        kernel( rankOffset, lengthTolerance, elemDofKey, faceDofKey, nodeManager, faceManager,
                subRegion, fluid, permeability, dt, localMatrix, localRhs );
        ElementBasedAssemblyKernel< NUM_FACES, IP >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    } );
  }

};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 * @brief Computes the norm of the face-based constitutive residuals, normalized by a
 *        pressure scale built from the adjacent cell pressures at the previous converged step.
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure_n >;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      SortedArrayView< localIndex const > const & regionFilter,
                      FaceManager const & faceManager,
                      SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_regionFilter( regionFilter ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_pres_n( singlePhaseFlowAccessors.get( fields::flow::pressure_n {} ) )
  {}

  GEOS_HOST_DEVICE
  void computePressureNormalizer( localIndex const kf,
                                  real64 & pressureNormalizer ) const
  {
    integer elemCounter = 0;
    for( integer k = 0; k < m_elemRegionList.size( 1 ); ++k )
    {
      localIndex const er  = m_elemRegionList[kf][k];
      localIndex const esr = m_elemSubRegionList[kf][k];
      localIndex const ei  = m_elemList[kf][k];
      bool const onBoundary = (er == -1 || esr == -1 || ei == -1);

      if( !onBoundary && m_regionFilter.contains( er ) )
      {
        pressureNormalizer = pressureNormalizer + LvArray::math::abs( m_pres_n[er][esr][ei] );
        elemCounter++;
      }
    }
    pressureNormalizer /= elemCounter;
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const kf,
                            LinfStackVariables & stack ) const override
  {
    if( m_dofNumber[kf] >= 0 )
    {
      real64 pressureNormalizer = 0.0;
      computePressureNormalizer( kf, pressureNormalizer );

      stack.localValue[0] = stack.localValue[0]
                            + LvArray::math::abs( m_localResidual[stack.localRow] ) / LvArray::math::max( m_minNormalizer, pressureNormalizer );
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const kf,
                          L2StackVariables & stack ) const override
  {
    if( m_dofNumber[kf] >= 0 )
    {
      real64 pressureNormalizer = 0.0;
      computePressureNormalizer( kf, pressureNormalizer );

      real64 const val = m_localResidual[stack.localRow];
      stack.localValue[0] += val * val;
      stack.localNormalizer[0] += LvArray::math::max( m_minNormalizer, pressureNormalizer );
    }
  }

protected:

  /// Filter to identify the target regions of the solver
  SortedArrayView< localIndex const > const m_regionFilter;

  /// Views on the maps face to elements
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// View on pressure at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   SortedArrayView< localIndex const > const & regionFilter,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   FaceManager const & faceManager,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = faceManager.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = faceManager.ghostRank();

    using kernelType = ResidualNormKernel;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               regionFilter, faceManager, flowAccessors, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( faceManager.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( faceManager.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

} // namespace singlePhaseMixedMFDKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP
