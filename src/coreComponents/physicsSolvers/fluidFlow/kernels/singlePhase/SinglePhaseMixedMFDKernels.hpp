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
 *
 * Flow equations of the mixed mimetic single-phase solver. Unknowns: face mass flux m_f and cell pressure p_K.
 * With sigma_f the orientation of the face f relative to the cell K, M the inner product weighted by K^{-1},
 * gamma = g . x the gravity coefficient and pi_f the pressure trace:
 *
 *   Darcy law      ( mu / rho )_K sum_j M_ij sigma_j m_j - p_K + rho_K ( gamma_K - gamma_f ) + pi_f = 0
 *   mass balance   ( phi rho V )_K - ( phi rho V )_K^n + dt sum_f sigma_f m_f = 0
 *
 * M = chi M_mfd + ( 1 - chi ) M_tpfa, chi the classification of the cell. A face whose cells all have chi = 0 is
 * condensed: its Darcy law reduces to the two-point relation A m_f = Phi_L - Phi_R, Phi = p - rho ( gamma_K - gamma_f ).
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "finiteVolume/mimeticInnerProducts/AdaptiveInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "mixedMimetic/MixedMimeticDispatch.hpp"
#include "mixedMimetic/MixedMimeticBoundaryConditions.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/consistency/ConsistencyAdaptationKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{
namespace singlePhaseMixedMFDKernels
{

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_FACE number of faces per element
 * @tparam IP the consistent inner product M_mfd
 * @brief Cell contributions: Darcy law on the non-condensed faces and mass balance.
 *
 * The Darcy law of a face is multiplied by sigma_f: the contributions of its two cells add up to the face
 * equation, in which pi_f cancels.
 */
template< integer NUM_FACE, typename IP >
class ElementBasedAssemblyKernel
{
public:

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using LocalToGlobalAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

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
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_myElemLocalToGlobal( subRegion.localToGlobalMap() ),
    m_faceOrientationCell( faceManager.getField< fields::mixedMimetic::faceOrientationCell >() ),
    m_nodePosition( nodeManager.referencePosition() ),
    m_elemPerm( permeability.permeability() ),
    m_elemPres( subRegion.getField< fields::flow::pressure >() ),
    m_faceFlux( faceManager.getField< fields::mixedMimetic::faceMassFlux >() ),
    m_faceArea( faceManager.faceArea() ),
    m_flowBc{ faceManager.getField< fields::mixedMimetic::flowBoundaryType >(),
             faceManager.getField< fields::mixedMimetic::flowBoundaryValue >(),
             faceManager.getField< fields::mixedMimetic::flowBoundaryCoefficient >() },
    m_faceStencilLabel( faceManager.getField< fields::mixedMimetic::faceStencilLabel >() ),
    m_mfdFlag( subRegion.getField< fields::mixedMimetic::mfdFlag >() ),
    m_elemDens( fluid.density() ),
    m_dElemDens( fluid.dDensity() ),
    m_mob( subRegion.getField< fields::flow::mobility >() ),
    m_dMob( subRegion.getField< fields::flow::dMobility >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Stack variables of the flow equations
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
      : massMatrix( NUM_FACE, NUM_FACE )
    {}

    /// M = chi M_mfd + ( 1 - chi ) M_tpfa
    stackArray2d< real64, NUM_FACE *NUM_FACE > massMatrix;

    /// sigma_i
    real64 orientation[NUM_FACE]{};

    /// 1 if m_f is prescribed (Neumann)
    integer isEssentialFace[NUM_FACE]{};

    /// 1 if m_f is condensed into a two-point relation
    integer isCondensedFace[NUM_FACE]{};

    /// sigma_i m_i, d( sigma_i m_i )/dm_i (0 if prescribed), ( M sigma m )_i and gamma_K - gamma_f
    real64 localFlux[NUM_FACE]{};
    real64 dLocalFlux[NUM_FACE]{};
    real64 mDotFlux[NUM_FACE]{};
    real64 gravCoefDif[NUM_FACE]{};

    /// Darcy residuals and their derivatives
    real64 faceResidual[NUM_FACE]{};
    real64 dFaceResidual_dFlux[NUM_FACE][NUM_FACE]{};
    real64 dFaceResidual_dPres[NUM_FACE]{};

    /// dt sum_f sigma_f m_f and its derivatives
    real64 divFlux = 0.0;
    real64 dDivFlux_dFlux[NUM_FACE]{};

    localIndex cellCenteredEqnRowIndex = 0;
    localIndex faceCenteredEqnRowIndex[NUM_FACE]{};
    globalIndex elemDofColIndex = 0;
    globalIndex faceDofColIndices[NUM_FACE]{};
  };

  /**
   * @brief Row and column indices of the cell.
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
   * @brief Darcy residuals and flux divergence of the cell.
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    real64 const perm[ 3 ] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };

    // ( mu / rho )_K and its pressure derivative
    real64 const mob = m_mob[ei];
    real64 const dMob_dPres = m_dMob[ei][DerivOffset::dP];
    real64 const invMob = 1.0 / mob;
    real64 const dInvMob_dPres = -dMob_dPres * invMob * invMob;

    real64 const ccDens = m_elemDens[ei][0];
    real64 const dCcDens_dPres = m_dElemDens[ei][0][DerivOffset::dP];

    // sigma_f = +1 for the cell of smallest global index of the face, -1 for the other. That cell is stored by
    // the owner of the face and synchronized: sigma_f is the same on every rank
    real64 ( &localFlux )[NUM_FACE] = stack.localFlux;
    globalIndex const myGlobalElem = m_myElemLocalToGlobal[ei];
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];

      stack.orientation[i] = ( myGlobalElem == m_faceOrientationCell[kf] ) ? 1.0 : -1.0;
      stack.gravCoefDif[i] = m_elemGravCoef[ei] - m_faceGravCoef[kf];

      stack.isEssentialFace[i] = m_flowBc[kf].isEssential() ? 1 : 0;
      stack.isCondensedFace[i] = ( m_faceStencilLabel[kf] == 0 && stack.isEssentialFace[i] == 0 ) ? 1 : 0;

      // prescribed flux: sigma m_f = g |f|, d( sigma m_f )/dm_f = 0
      real64 const faceFlux = stack.isEssentialFace[i] == 1 ? m_flowBc[kf].essentialFlux( stack.orientation[i], m_faceArea[kf] )
                                                            : m_faceFlux[kf];
      localFlux[i] = stack.orientation[i] * faceFlux;
      stack.dLocalFlux[i] = stack.isEssentialFace[i] == 1 ? 0.0 : stack.orientation[i];
    }

    // M, needed if one face of the cell is non-condensed
    bool anyLiveFace = false;
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      anyLiveFace = anyLiveFace || ( stack.isEssentialFace[i] == 0 && stack.isCondensedFace[i] == 0 );
    }
    if( anyLiveFace )
    {
      real64 const chi = static_cast< real64 >( m_mfdFlag[ei] );
      mimeticInnerProduct::AdaptiveInnerProduct< IP >::template computeM< NUM_FACE >( m_nodePosition,
                                                                                      m_faceToNodes,
                                                                                      m_elemToFaces[ei],
                                                                                      m_elemCenter[ei],
                                                                                      m_elemVolume[ei],
                                                                                      perm,
                                                                                      m_lengthTolerance,
                                                                                      chi,
                                                                                      stack.massMatrix );
    }

    // Darcy law on the non-condensed faces: ( mu / rho )_K sum_j M_ij sigma_j m_j - p_K + rho_K ( gamma_K - gamma_f ) + pi_f = 0.
    // pi_f cancels between the two cells of an interior face; the boundary condition defines it on a boundary face
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      if( stack.isEssentialFace[i] == 1 || stack.isCondensedFace[i] == 1 )
      {
        continue;
      }

      localIndex const kf = m_elemToFaces[ei][i];

      real64 mDotFlux = 0.0;
      for( integer j = 0; j < NUM_FACE; ++j )
      {
        mDotFlux += stack.massMatrix( i, j ) * localFlux[j];
        stack.dFaceResidual_dFlux[i][j] = invMob * stack.massMatrix( i, j ) * stack.dLocalFlux[j];
      }
      stack.mDotFlux[i] = mDotFlux;

      real64 const gravCoefDif = stack.gravCoefDif[i];

      stack.faceResidual[i] = invMob * mDotFlux - m_elemPres[ei] + ccDens * gravCoefDif;
      stack.dFaceResidual_dPres[i] = dInvMob_dPres * mDotFlux - 1.0 + dCcDens_dPres * gravCoefDif;

      mixedMimeticBoundary::FaceBoundaryCondition const bc = m_flowBc[kf];
      if( bc.type != mixedMimeticBoundary::BoundaryType::interior )
      {
        real64 dTrace_dFlux = 0.0;
        stack.faceResidual[i] += bc.trace( stack.orientation[i], m_faceArea[kf], localFlux[i], dTrace_dFlux );
        stack.dFaceResidual_dFlux[i][i] += dTrace_dFlux * stack.orientation[i];
      }
    }

    // dt sum_f sigma_f m_f over the non-condensed faces
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      if( stack.isCondensedFace[i] == 1 )
      {
        stack.dDivFlux_dFlux[i] = 0.0;
        continue;
      }
      stack.divFlux += m_dt * localFlux[i];
      stack.dDivFlux_dFlux[i] = m_dt * stack.dLocalFlux[i];
    }
  }

  /**
   * @brief Add the contributions of the cell to the system.
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // mass balance of the cell
    if( m_elemGhostRank[ei] < 0 )
    {
      m_localRhs[stack.cellCenteredEqnRowIndex] += stack.divFlux;

      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( stack.cellCenteredEqnRowIndex,
                                                                  &stack.faceDofColIndices[0],
                                                                  &stack.dDivFlux_dFlux[0],
                                                                  NUM_FACE );
    }

    // Darcy law, multiplied by sigma_f
    for( integer i = 0; i < NUM_FACE; ++i )
    {
      localIndex const kf = m_elemToFaces[ei][i];
      if( m_faceGhostRank[kf] >= 0 || stack.isCondensedFace[i] == 1 )
      {
        continue;
      }

      if( stack.isEssentialFace[i] == 1 )
      {
        // Neumann condition: m_f - sigma g |f| = 0
        real64 const one = 1.0;
        real64 const target = m_flowBc[kf].essentialFlux( stack.orientation[i], m_faceArea[kf] );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[stack.faceCenteredEqnRowIndex[i]], m_faceFlux[kf] - target );
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
   * @brief Launch the kernel over the cells.
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

  /// first global row of the rank
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
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;
  arrayView1d< globalIndex const > const m_myElemLocalToGlobal;
  /// global index of the cell with sigma_f = +1
  arrayView1d< globalIndex const > const m_faceOrientationCell;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;

  /// permeability
  arrayView3d< real64 const > const m_elemPerm;

  /// primary variables and boundary data
  arrayView1d< real64 const > const m_elemPres;
  arrayView1d< real64 const > const m_faceFlux;
  arrayView1d< real64 const > const m_faceArea;
  mixedMimeticBoundary::FaceBoundaryView const m_flowBc;
  arrayView1d< integer const > const m_faceStencilLabel;

  /// chi
  arrayView1d< integer const > const m_mfdFlag;

  /// fluid data
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_elemDens;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dElemDens;
  arrayView1d< real64 const > const m_mob;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_dMob;

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

/******************************** TpfaCondensedFluxKernel ********************************/

/**
 * @class TpfaCondensedFluxKernel
 * @brief Contributions of the condensed faces, with L and R the cells of the face. M_tpfa is diagonal, so the
 *        Darcy laws of L and R give m_f = ( Phi_L - Phi_R ) / A, A = sum_K ( mu / rho )_K / t_K, t_K the one-sided
 *        transmissibility. The kernel adds +/- dt m_f( p_L, p_R ) to the mass balances of L and R and the face
 *        equation A m_f = Phi_L - Phi_R. No other equation depends on m_f: the system is equivalent to the
 *        non-condensed one, and m_f is decoupled from the cell unknowns.
 */
class TpfaCondensedFluxKernel
{
public:

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using FlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::flow::mobility,
                      fields::flow::dMobility >;

  using FluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using GhostRankAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > >;

  TpfaCondensedFluxKernel( globalIndex const rankOffset,
                           real64 const & lengthTolerance,
                           string const faceDofKey,
                           NodeManager const & nodeManager,
                           FaceManager const & faceManager,
                           DofNumberAccessor const & elemDofNumber,
                           DofNumberAccessor const & elemLocalToGlobal,
                           GhostRankAccessor const & elemGhostRank,
                           ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const & elemCenter,
                           FlowAccessors const & flowAccessors,
                           FluidAccessors const & fluidAccessors,
                           PermeabilityAccessors const & permAccessors,
                           SortedArrayView< localIndex const > const & regionFilter,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_lengthTolerance( lengthTolerance ),
    m_dt( dt ),
    m_faceGhostRank( faceManager.ghostRank() ),
    m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
    m_faceToNodes( faceManager.nodeList().toViewConst() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_faceFlux( faceManager.getField< fields::mixedMimetic::faceMassFlux >() ),
    m_faceArea( faceManager.faceArea() ),
    m_flowBc{ faceManager.getField< fields::mixedMimetic::flowBoundaryType >(),
             faceManager.getField< fields::mixedMimetic::flowBoundaryValue >(),
             faceManager.getField< fields::mixedMimetic::flowBoundaryCoefficient >() },
    m_faceStencilLabel( faceManager.getField< fields::mixedMimetic::faceStencilLabel >() ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_nodePosition( nodeManager.referencePosition() ),
    m_regionFilter( regionFilter ),
    m_elemDofNumber( elemDofNumber.toNestedViewConst() ),
    m_elemLocalToGlobal( elemLocalToGlobal.toNestedViewConst() ),
    m_elemGhostRank( elemGhostRank.toNestedViewConst() ),
    m_elemCenter( elemCenter.toNestedViewConst() ),
    m_pres( flowAccessors.get( fields::flow::pressure {} ) ),
    m_elemGravCoef( flowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_mob( flowAccessors.get( fields::flow::mobility {} ) ),
    m_dMob( flowAccessors.get( fields::flow::dMobility {} ) ),
    m_dens( fluidAccessors.get( fields::singlefluid::density {} ) ),
    m_dDens( fluidAccessors.get( fields::singlefluid::dDensity {} ) ),
    m_elemPerm( permAccessors.get( fields::permeability::permeability {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @brief Cells of a face in the target regions, the cell with sigma_f = +1 first.
   * @return the number of cells
   */
  GEOS_HOST_DEVICE
  integer gatherCells( localIndex const kf,
                       localIndex ( & er )[2],
                       localIndex ( & esr )[2],
                       localIndex ( & ei )[2] ) const
  {
    integer numElems = 0;
    for( integer k = 0; k < m_elemRegionList.size( 1 ); ++k )
    {
      localIndex const erk  = m_elemRegionList[kf][k];
      localIndex const esrk = m_elemSubRegionList[kf][k];
      localIndex const eik  = m_elemList[kf][k];
      if( erk >= 0 && esrk >= 0 && eik >= 0 && m_regionFilter.contains( erk ) )
      {
        er[numElems] = erk; esr[numElems] = esrk; ei[numElems] = eik;
        numElems++;
      }
    }
    // L is the cell of smallest global index: sigma_f = +1, as in the cell kernel
    if( numElems == 2 &&
        m_elemLocalToGlobal[er[1]][esr[1]][ei[1]] < m_elemLocalToGlobal[er[0]][esr[0]][ei[0]] )
    {
      localIndex tmp;
      tmp = er[0]; er[0] = er[1]; er[1] = tmp;
      tmp = esr[0]; esr[0] = esr[1]; esr[1] = tmp;
      tmp = ei[0]; ei[0] = ei[1]; ei[1] = tmp;
    }
    return numElems;
  }

  /// A, pot = Phi_L - Phi_R, m_f = pot / A and their derivatives with respect to the cell unknowns
  template< integer NUM_DERIV >
  struct TwoPointFlux
  {
    real64 A = 0.0;
    real64 pot = 0.0;
    real64 flux = 0.0;
    real64 dA[2][NUM_DERIV]{};
    real64 dPot[2][NUM_DERIV]{};
    real64 dFlux[2][NUM_DERIV]{};
  };

  /**
   * @brief One-sided transmissibility t_K of a cell through a face for a diagonal coefficient, as in TPFAInnerProduct::computeM.
   */
  GEOS_HOST_DEVICE
  real64 oneSidedTrans( localIndex const kf,
                        localIndex const er, localIndex const esr, localIndex const ei,
                        real64 const (&coef)[3] ) const
  {
    real64 const areaTolerance = m_lengthTolerance * m_lengthTolerance;
    real64 const weightTolerance = 1e-30 * m_lengthTolerance;
    return LvArray::math::max( mimeticInnerProduct::TPFAInnerProduct::computeOneSidedTrans( m_nodePosition,
                                                                                            m_faceToNodes,
                                                                                            kf,
                                                                                            m_elemCenter[er][esr][ei],
                                                                                            coef,
                                                                                            areaTolerance ),
                               weightTolerance );
  }

  /**
   * @brief Two-point mass flux m_f = ( Phi_L - Phi_R ) / A, A = sum_K ( mu / rho )_K / t_K. On a boundary face
   *        Phi_R is the boundary value and A includes 1 / ( alpha |f| ).
   * @tparam NUM_DERIV 1 (pressure) or 2 (pressure and temperature) derivatives of the cell unknowns
   */
  template< integer NUM_DERIV >
  GEOS_HOST_DEVICE
  void computeTwoPointMassFlux( localIndex const kf,
                                localIndex const (&er)[2], localIndex const (&esr)[2], localIndex const (&ei)[2],
                                integer const numElems,
                                TwoPointFlux< NUM_DERIV > & tp ) const
  {
    real64 grav[2]{}, dGrav[2][NUM_DERIV]{};
    for( integer k = 0; k < numElems; ++k )
    {
      real64 const perm[3] = { m_elemPerm[er[k]][esr[k]][ei[k]][0][0],
                               m_elemPerm[er[k]][esr[k]][ei[k]][0][1],
                               m_elemPerm[er[k]][esr[k]][ei[k]][0][2] };
      real64 const t = oneSidedTrans( kf, er[k], esr[k], ei[k], perm );
      real64 const invMob = 1.0 / m_mob[er[k]][esr[k]][ei[k]];
      real64 const gravCoefDif = m_elemGravCoef[er[k]][esr[k]][ei[k]] - m_faceGravCoef[kf];
      tp.A += invMob / t;
      grav[k] = m_dens[er[k]][esr[k]][ei[k]][0] * gravCoefDif;
      for( integer d = 0; d < NUM_DERIV; ++d )
      {
        tp.dA[k][d] = -m_dMob[er[k]][esr[k]][ei[k]][d] * invMob * invMob / t;
        dGrav[k][d] = m_dDens[er[k]][esr[k]][ei[k]][0][d] * gravCoefDif;
      }
    }
    if( numElems == 1 )
    {
      tp.A += m_flowBc[kf].resistance( m_faceArea[kf] );
    }
    real64 const pL = m_pres[er[0]][esr[0]][ei[0]];
    real64 const pR = ( numElems == 2 ) ? m_pres[er[1]][esr[1]][ei[1]] : m_flowBc[kf].value;
    tp.pot = ( pL - pR ) - grav[0] + ( ( numElems == 2 ) ? grav[1] : 0.0 );
    tp.flux = tp.pot / tp.A;
    for( integer k = 0; k < numElems; ++k )
    {
      real64 const sign = ( k == 0 ) ? 1.0 : -1.0;
      for( integer d = 0; d < NUM_DERIV; ++d )
      {
        tp.dPot[k][d] = ( d == 0 ? sign : 0.0 ) - sign * dGrav[k][d];
        tp.dFlux[k][d] = ( tp.dPot[k][d] - tp.flux * tp.dA[k][d] ) / tp.A;
      }
    }
  }

  GEOS_HOST_DEVICE
  void compute( localIndex const kf ) const
  {
    if( m_faceStencilLabel[kf] != 0 )
    {
      return;
    }

    localIndex er[2]{}, esr[2]{}, ei[2]{};
    integer const numElems = gatherCells( kf, er, esr, ei );
    if( numElems == 0 )
    {
      return;
    }

    if( m_flowBc[kf].isEssential() )
    {
      // prescribed flux: assembled by the cell kernel
      return;
    }

    TwoPointFlux< 1 > tp;
    computeTwoPointMassFlux< 1 >( kf, er, esr, ei, numElems, tp );
    real64 const A = tp.A;
    real64 const pot = tp.pot;
    real64 const mFlux = tp.flux;
    real64 const dA_dp[2] = { tp.dA[0][0], tp.dA[1][0] };
    real64 const dPot_dp[2] = { tp.dPot[0][0], tp.dPot[1][0] };
    real64 const dFlux_dp[2] = { tp.dFlux[0][0], tp.dFlux[1][0] };

    globalIndex elemDofCols[2];
    for( integer k = 0; k < numElems; ++k )
    {
      elemDofCols[k] = m_elemDofNumber[er[k]][esr[k]][ei[k]];
    }

    // mass balances of L and R: +/- dt m_f
    for( integer k = 0; k < numElems; ++k )
    {
      if( m_elemGhostRank[er[k]][esr[k]][ei[k]] >= 0 )
      {
        continue;
      }
      real64 const sign = ( k == 0 ) ? 1.0 : -1.0;
      localIndex const row = m_elemDofNumber[er[k]][esr[k]][ei[k]] - m_rankOffset;

      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[row], sign * m_dt * mFlux );

      real64 vals[2];
      for( integer j = 0; j < numElems; ++j )
      {
        vals[j] = sign * m_dt * dFlux_dp[j];
      }
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( row,
                                                                          elemDofCols,
                                                                          vals,
                                                                          numElems );
    }

    // face equation: A m_f - ( Phi_L - Phi_R ) = 0
    if( m_faceGhostRank[kf] < 0 )
    {
      localIndex const row = m_faceDofNumber[kf] - m_rankOffset;
      globalIndex const faceDofCol = m_faceDofNumber[kf];
      real64 const mCurrent = m_faceFlux[kf];

      m_localRhs[row] += A * mCurrent - pot;

      m_localMatrix.addToRow< serialAtomic >( row, &faceDofCol, &A, 1 );

      real64 vals[2];
      for( integer k = 0; k < numElems; ++k )
      {
        vals[k] = dA_dp[k] * mCurrent - dPot_dp[k];
      }
      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( row,
                                                                  elemDofCols,
                                                                  vals,
                                                                  numElems );
    }
  }

  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numFaces,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numFaces, [=] GEOS_HOST_DEVICE ( localIndex const kf )
    {
      kernelComponent.compute( kf );
    } );
  }

protected:

  globalIndex const m_rankOffset;
  real64 const m_lengthTolerance;
  real64 const m_dt;

  /// face data
  arrayView1d< integer const > const m_faceGhostRank;
  arrayView1d< globalIndex const > const m_faceDofNumber;
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  arrayView1d< real64 const > const m_faceGravCoef;
  arrayView1d< real64 const > const m_faceFlux;
  arrayView1d< real64 const > const m_faceArea;
  mixedMimeticBoundary::FaceBoundaryView const m_flowBc;
  arrayView1d< integer const > const m_faceStencilLabel;
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;
  SortedArrayView< localIndex const > const m_regionFilter;

  /// cell data
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumber;
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemLocalToGlobal;
  ElementViewConst< arrayView1d< integer const > > const m_elemGhostRank;
  ElementViewConst< arrayView2d< real64 const > > const m_elemCenter;
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView1d< real64 const > > const m_elemGravCoef;
  ElementViewConst< arrayView1d< real64 const > > const m_mob;
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_dMob;
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_dens;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dDens;
  ElementViewConst< arrayView3d< real64 const > > const m_elemPerm;

  /// system views
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  arrayView1d< real64 > const m_localRhs;

};

/**
 * @class TpfaCondensedFluxKernelFactory
 */
class TpfaCondensedFluxKernelFactory
{
public:

  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   string const solverName,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   SortedArrayView< localIndex const > const & regionFilter,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    TpfaCondensedFluxKernel::DofNumberAccessor elemDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( solverName + "/accessors/" + elemDofKey );

    TpfaCondensedFluxKernel::DofNumberAccessor const elemLocalToGlobal =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );

    TpfaCondensedFluxKernel::GhostRankAccessor const elemGhostRank =
      elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
      elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

    TpfaCondensedFluxKernel::FlowAccessors flowAccessors( elemManager, solverName );
    TpfaCondensedFluxKernel::FluidAccessors fluidAccessors( elemManager, solverName );
    TpfaCondensedFluxKernel::PermeabilityAccessors permAccessors( elemManager, solverName );

    TpfaCondensedFluxKernel kernel( rankOffset, lengthTolerance, faceDofKey,
                                    nodeManager, faceManager,
                                    elemDofNumber, elemLocalToGlobal, elemGhostRank, elemCenter,
                                    flowAccessors, fluidAccessors, permAccessors,
                                    regionFilter, dt, localMatrix, localRhs );

    TpfaCondensedFluxKernel::launch< POLICY >( faceManager.size(), kernel );
  }

};

} // namespace singlePhaseMixedMFDKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MIXEDMFDKERNELS_HPP
