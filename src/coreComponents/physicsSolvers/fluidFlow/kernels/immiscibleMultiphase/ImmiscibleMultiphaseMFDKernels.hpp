/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * (c) GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ImmiscibleMultiphaseMFDKernels.hpp
 *
 * Hybrid Mimetic (mixed) finite volume kernels for immiscible multiphase flow pressure system.
 * This initial implementation mirrors the single-phase HybridFVM structure and reuses
 * a summed (upwind) total mobility with an effective density (currently phase 0 only).
 *
 * TODO: extend to proper phase-weighted density and face-based upwinding per phase.
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASE_MFDKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASE_MFDKERNELS_HPP

#include "common/DataTypes.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
// +++ added for accumulation (porosity access)
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "mesh/MeshLevel.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/kernels/HybridFVMHelperKernels.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"

namespace geos
{
namespace immiscibleMultiphaseMFDKernels
{

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  switch( value )
  {
    case 4:  lambda( std::integral_constant< T, 4 >() ); return;
    case 5:  lambda( std::integral_constant< T, 5 >() ); return;
    case 6:  lambda( std::integral_constant< T, 6 >() ); return;
    case 7:  lambda( std::integral_constant< T, 7 >() ); return;
    case 8:  lambda( std::integral_constant< T, 8 >() ); return;
    case 9:  lambda( std::integral_constant< T, 9 >() ); return;
    case 10: lambda( std::integral_constant< T, 10 >() ); return;
    case 11: lambda( std::integral_constant< T, 11 >() ); return;
    case 12: lambda( std::integral_constant< T, 12 >() ); return;
    case 13: lambda( std::integral_constant< T, 13 >() ); return;
    default: GEOS_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

/**
 * AveragePressureGradientKernel (reuse for diagnostics)
 */
template< integer NUM_FACES >
class AveragePressureGradientKernel
{
public:
  AveragePressureGradientKernel( CellElementSubRegion & subRegion,
                                 FaceManager const & faceManager )
    : m_facePressure( faceManager.getField< fields::flow::facePressure >() ),
      m_faceCenter( faceManager.faceCenter() ),
      m_pres( subRegion.getField< fields::flow::pressure >() ),
      m_elemCenter( subRegion.getElementCenter() ),
      m_elemsToFaces( subRegion.faceList() ),
      m_presGradient( subRegion.getField< fields::flow::pressureGradient >() )
  {}

  struct StackVariables
  {
    GEOS_HOST_DEVICE StackVariables(): coordinates( NUM_FACES+1, 4 ), pressures( NUM_FACES+1 ), presGradientLocal( 4 ) {}
    stackArray2d< real64, (NUM_FACES + 1) * 4 > coordinates;
    stackArray1d< real64, NUM_FACES + 1 > pressures;
    stackArray1d< real64, 4 > presGradientLocal;
  };

  GEOS_HOST_DEVICE
  void compute( localIndex const ei, StackVariables & stack ) const
  {
    for( integer d=0; d<3; ++d )
    {
      stack.coordinates( 0, d ) = m_elemCenter( ei, d );
    }
    stack.coordinates( 0, 3 ) = 1.0;
    stack.pressures[0] = m_pres[ei];
    for( integer f=0; f<NUM_FACES; ++f )
    {
      localIndex const lf = m_elemsToFaces( ei, f );
      stack.pressures[f+1] = m_facePressure[lf];
      for( integer d=0; d<3; ++d ) stack.coordinates( f+1, d ) = m_faceCenter[lf][d];
      stack.coordinates( f+1, 3 ) = 1.0;
    }
    BlasLapackLA::matrixLeastSquaresSolutionSolve( stack.coordinates, stack.pressures, stack.presGradientLocal );
    for( integer d=0; d<3; ++d ) m_presGradient( ei, d ) = stack.presGradientLocal[d];
  }

  template< typename POLICY, typename KERNEL_TYPE >
  static void launch( localIndex const n, KERNEL_TYPE const & kern )
  {
    forAll< POLICY >( n, [=] ( localIndex const ei )
    {
      typename KERNEL_TYPE::StackVariables stack; kern.compute( ei, stack );
    } );
  }
private:
  arrayView1d< real64 const > const m_facePressure;
  arrayView2d< real64 const > const m_faceCenter;
  arrayView1d< real64 const > const m_pres;
  arrayView2d< real64 const > const m_elemCenter;
  arrayView2d< localIndex const > const m_elemsToFaces;
  arrayView2d< real64 > const m_presGradient;
};

class AveragePressureGradientKernelFactory
{
public:
  template< typename POLICY >
  static void createAndLaunch( CellElementSubRegion & subRegion, FaceManager const & faceManager )
  {
    internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
    {
      AveragePressureGradientKernel< NF > k( subRegion, faceManager );
      AveragePressureGradientKernel< NF >::template launch< POLICY >( subRegion.size(), k );
    } );
  }
};

/**
 * Element-based assembly kernel for hybrid pressure system (immiscible multiphase)
 * Currently uses an effective density = density of phase 0.
 */
template< integer NUM_FACE, typename IP >
class ElementBasedAssemblyKernel
{
public:
  using DerivMob = immiscibleFlow::DerivativeOffset; // for mobility derivatives
  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;
  using SatAccessor = ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >;
  using PhaseMobAccessor = ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > >;
  using DPhaseMobAccessor = ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, immiscibleFlow::USD_PHASE_DS > >;
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              localIndex const er,
                              localIndex const esr,
                              real64 const & lengthTolerance,
                              string const faceDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              CellElementSubRegion const & subRegion,
                              DofNumberAccessor const & elemDofNumberAccessor,
                              SatAccessor const & phaseVolFracAccessor,
                              constitutive::TwoPhaseImmiscibleFluid const & fluid,
                              constitutive::PermeabilityBase const & permeability,
                              SortedArrayView< localIndex const > const & regionFilter,
                              integer const indepPhaseIndex,
                              real64 const & dt,
                              bool const assembleCellEq,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              PhaseMobAccessor const & phaseMobAccessor,
                              DPhaseMobAccessor const & dPhaseMobAccessor )
    : m_rankOffset( rankOffset ), m_er( er ), m_esr( esr ), m_lengthTolerance( lengthTolerance ), m_dt( dt ),
      m_elemGhostRank( subRegion.ghostRank() ),
      m_elemDofNumber( elemDofNumberAccessor.toNestedViewConst() ),
      m_phaseVolFracAll( phaseVolFracAccessor.toNestedViewConst() ),
      m_faceGhostRank( faceManager.ghostRank() ), m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
      m_elemToFaces( subRegion.faceList().toViewConst() ), m_elemCenter( subRegion.getElementCenter() ), m_elemVolume( subRegion.getElementVolume() ), m_elemGravCoef( subRegion.getField< fields::flow::gravityCoefficient >() ),
      m_faceToNodes( faceManager.nodeList().toViewConst() ), m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ), m_regionFilter( regionFilter ),
      m_nodePosition( nodeManager.referencePosition() ), m_elemRegionList( faceManager.elementRegionList() ), m_elemSubRegionList( faceManager.elementSubRegionList() ), m_elemList( faceManager.elementList() ),
      m_elemPerm( permeability.permeability() ), m_transMultiplier( faceManager.getField< fields::flow::transMultiplier >() ),
      m_elemPres( subRegion.getField< fields::flow::pressure >() ), m_facePres( faceManager.getField< fields::flow::facePressure >() ),
      m_phaseDens( fluid.phaseDensity() ), m_dPhaseDens( fluid.dPhaseDensity() ),
      m_phaseMobAll( phaseMobAccessor.toNestedViewConst() ), m_dPhaseMobAll( dPhaseMobAccessor.toNestedViewConst() ),
      m_assembleCellEquation( assembleCellEq ), m_indep( indepPhaseIndex ),
      m_localMatrix( localMatrix ), m_localRhs( localRhs )
  {}

  struct StackVariables
  {
    GEOS_HOST_DEVICE StackVariables(): transMatrix( NUM_FACE, NUM_FACE ) {}
    stackArray2d< real64, NUM_FACE * NUM_FACE > transMatrix;
    // Mass fluxes (one-sided per face i)
    real64 MassFlux[NUM_FACE]{};
    real64 dMassFlux_dPres[NUM_FACE]{};
    real64 dMassFlux_dS[NUM_FACE]{};
    real64 dMassFlux_dFacePres[NUM_FACE][NUM_FACE]{};
    // Pressure equation: divergence of total mass flux
    real64 divMassFluxes = 0;
    real64 dDivMassFluxes_dP = 0.0;
    real64 dDivMassFluxes_dS = 0.0;
    real64 dDivMassFluxes_dFaceVars[NUM_FACE]{};
    // Saturation equation: divergence of upwinded transported saturation with mass flux
    real64 divSatFluxes = 0.0;
    real64 dDivSatFluxes_dP = 0.0;
    real64 dDivSatFluxes_dS = 0.0;
    real64 dDivSatFluxes_dFaceVars[NUM_FACE]{};
    // Neighbor saturation couplings (via upwind when inflow)
    localIndex numNeiCols = 0;
    globalIndex neiCols[2*NUM_FACE]{}; // allow both P and S per neighbor face
    real64 neiVals[2*NUM_FACE]{};
    // Row/col bookkeeping
    localIndex cellRow = 0;
    localIndex faceRow[NUM_FACE]{};
    globalIndex elemCols[2]{}; // [P, S_indep]
    globalIndex faceCols[NUM_FACE]{};
  };

  GEOS_HOST_DEVICE void setup( localIndex const ei, StackVariables & s ) const
  {
    s.cellRow = m_elemDofNumber[m_er][m_esr][ei] - m_rankOffset;
    // Element columns: pressure and saturation of independent phase are consecutive
    s.elemCols[0] = m_elemDofNumber[m_er][m_esr][ei];
    s.elemCols[1] = s.elemCols[0] + 1;
    for( integer f=0; f<NUM_FACE; ++f )
    {
      localIndex const lf = m_elemToFaces[ei][f];
      s.faceRow[f] = m_faceDofNumber[lf] - m_rankOffset;
      s.faceCols[f] = m_faceDofNumber[lf];
    }
  }
  
  GEOS_HOST_DEVICE void computeOverallMassFlux( localIndex const ei, StackVariables & s ) const
  {
    using Deriv = DerivMob; // alias for readability
    for( integer i=0; i<NUM_FACE; ++i )
    {
      for( integer j=0; j<NUM_FACE; ++j )
      {
        // Local pressure and gravity terms (cell-centered and face values associated to j)
        real64 const ccPres = m_elemPres[ei];
        real64 const fPres = m_facePres[m_elemToFaces[ei][j]];
        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][j]];

        // Total mobility and its derivatives at element ei
        real64 Lambda = 0.0;
        real64 dLambda_dP = 0.0;
        real64 dLambda_dS = 0.0; // derivative w.r.t. independent saturation (apply sign map below)

        // Mixture density rho_mix = (sum_i rho_i * lambda_i) / (sum_i lambda_i)
        // Keep pressure and saturation derivatives since lambda depends on both
        real64 sumLambda = 0.0;
        real64 dSumLambda_dP = 0.0;
        real64 dSumLambda_dS = 0.0;
        real64 sumRhoLambda = 0.0;
        real64 dSumRhoLambda_dP = 0.0;
        real64 dSumRhoLambda_dS = 0.0;

        // two-phase for now
        for( integer ip=0; ip<2; ++ip )
        {
          real64 const rho = m_phaseDens[ei][0][ip];
          real64 const drho_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
          real64 const lambda = m_phaseMobAll[m_er][m_esr][ei][ip];
          real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dP];
          real64 const dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dS];

          // accumulate mobility and its derivs
          Lambda += lambda;
          dLambda_dP += dlambda_dP;
          dLambda_dS += dlambda_dS_raw; // adjust sign after loop if indep != 0

          // terms for mixture density and its derivs
          sumLambda += lambda;
          dSumLambda_dP += dlambda_dP;
          dSumLambda_dS += dlambda_dS_raw; // adjust sign after loop if indep != 0

          sumRhoLambda += rho * lambda;
          dSumRhoLambda_dP += drho_dP * lambda + rho * dlambda_dP;
          dSumRhoLambda_dS += /* drho/dS ~ 0 */ rho * dlambda_dS_raw; // adjust sign after loop if indep != 0
        }

        // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
        real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
        dLambda_dS *= sgnS;
        dSumLambda_dS *= sgnS;
        dSumRhoLambda_dS *= sgnS;

        // Safe denominator (total mobility is expected > 0, but guard for robustness)
        real64 const eps = 1e-30;
        real64 const denom = (fabs( sumLambda ) > eps) ? sumLambda : (sumLambda >= 0.0 ? eps : -eps);
        real64 const densMix = sumRhoLambda / denom;
        real64 const dDensMix_dP = ( dSumRhoLambda_dP * denom - sumRhoLambda * dSumLambda_dP ) / ( denom * denom );
        real64 const dDensMix_dS = ( dSumRhoLambda_dS * denom - sumRhoLambda * dSumLambda_dS ) / ( denom * denom );

        // Potential difference and its derivatives
        real64 const presDif = ccPres - fPres;
        real64 const gravCoefDif = ccGravCoef - fGravCoef;
        real64 const gravTerm = densMix * gravCoefDif;
        real64 const dGravTerm_dP = dDensMix_dP * gravCoefDif;
        real64 const dGravTerm_dS = dDensMix_dS * gravCoefDif;
        real64 const potDif = presDif - gravTerm;
        real64 const dPotDif_dP = 1.0 - dGravTerm_dP;
        real64 const dPotDif_dS = - dGravTerm_dS;
        real64 const dPotDif_dFaceP = -1.0;

        // Overall mass flux and derivatives
        real64 const T_ij = s.transMatrix[i][j];
        s.MassFlux[i] += Lambda * T_ij * potDif;
        s.dMassFlux_dPres[i] += T_ij * ( Lambda * dPotDif_dP + dLambda_dP * potDif );
        s.dMassFlux_dS[i] += T_ij * ( Lambda * dPotDif_dS + dLambda_dS * potDif );
        s.dMassFlux_dFacePres[i][j] += T_ij * ( Lambda * dPotDif_dFaceP ); // = - T_ij * Lambda
      }
    }
  }

  
  GEOS_HOST_DEVICE
  void computeOverallMassFluxDivergence( localIndex const ei, StackVariables & s ) const
  {
    for( integer i=0; i<NUM_FACE; ++i )
    {
      real64 const F = s.MassFlux[i];
      real64 const dF_dP = s.dMassFlux_dPres[i];
      real64 const dF_dS = s.dMassFlux_dS[i];
      // residual
      s.divMassFluxes += m_dt * F;
      // jacobians wrt element DOFs
      s.dDivMassFluxes_dP += m_dt * dF_dP;
      s.dDivMassFluxes_dS += m_dt * dF_dS;
      // wrt face pressures
      for( integer j=0; j<NUM_FACE; ++j ) s.dDivMassFluxes_dFaceVars[j] += m_dt * s.dMassFlux_dFacePres[i][j];
    }
  }

  GEOS_HOST_DEVICE
  void computeSaturationTransport( localIndex const ei, StackVariables & s ) const
  {
    // local fractional flow components for independent phase
    auto fracFlow = [this] GEOS_HOST_DEVICE ( localIndex const er,
                                              localIndex const esr,
                                              localIndex const ei_local,
                                              integer const indep,
                                              real64 & f,
                                              real64 & df_dP,
                                              real64 & df_dS )
    {
      using Deriv = DerivMob;
      // two-phase: independent ip = indep, dependent ip = 1 - indep
      integer const dep = 1 - indep;
      real64 const lam_ind = m_phaseMobAll[er][esr][ei_local][indep];
      real64 const lam_dep = m_phaseMobAll[er][esr][ei_local][dep];
      real64 const dlam_ind_dP = m_dPhaseMobAll[er][esr][ei_local][indep][Deriv::dP];
      real64 const dlam_dep_dP = m_dPhaseMobAll[er][esr][ei_local][dep][Deriv::dP];
      // Stored d/dS is w.r.t phase-0 saturation; map to configured independent saturation
      real64 const sgnS = ( indep == 0 ? 1.0 : -1.0 );
      real64 const dlam_ind_dS = sgnS * m_dPhaseMobAll[er][esr][ei_local][indep][Deriv::dS];
      real64 const dlam_dep_dS = sgnS * m_dPhaseMobAll[er][esr][ei_local][dep][Deriv::dS];

      real64 const Lambda = lam_ind + lam_dep;

      f = lam_ind / Lambda;
      // df/dx = (dlam_ind*Lambda - lam_ind*dLambda) / Lambda^2
      real64 const dLambda_dP = dlam_ind_dP + dlam_dep_dP;
      real64 const dLambda_dS = dlam_ind_dS + dlam_dep_dS;
      df_dP = ( dlam_ind_dP * Lambda - lam_ind * dLambda_dP ) / ( Lambda * Lambda );
      df_dS = ( dlam_ind_dS * Lambda - lam_ind * dLambda_dS ) / ( Lambda * Lambda );
      df_dS = sgnS * 1.0;
    };

    // Precompute local f and derivatives
    real64 f_loc = 0.0, df_loc_dP = 0.0, df_loc_dS = 0.0;
    fracFlow( m_er, m_esr, ei, m_indep, f_loc, df_loc_dP, df_loc_dS );

    for( integer i=0; i<NUM_FACE; ++i )
    {
      // Determine neighbor across face i
      localIndex const lf = m_elemToFaces[ei][i];
      // face-adjacent elements
      localIndex const er0 = m_elemRegionList[lf][0];
      localIndex const esr0 = m_elemSubRegionList[lf][0];
      localIndex const ei0 = m_elemList[lf][0];
      localIndex const er1 = m_elemRegionList[lf][1];
      localIndex const esr1 = m_elemSubRegionList[lf][1];
      localIndex const ei1 = m_elemList[lf][1];

      // identify neighbor indices
      bool const side0IsSelf = (er0 == m_er && esr0 == m_esr && ei0 == ei);
      bool const side1IsSelf = (er1 == m_er && esr1 == m_esr && ei1 == ei);
      localIndex ner = -1, nesr = -1, nei = -1;
      if( side0IsSelf && er1 >= 0 ) { ner = er1; nesr = esr1; nei = ei1; }
      else if( side1IsSelf && er0 >= 0 ) { ner = er0; nesr = esr0; nei = ei0; }

      // mass flux and derivatives for this one-sided face i (from cell to face)
      real64 const F = s.MassFlux[i];
      real64 const dF_dP = s.dMassFlux_dPres[i];
      real64 const dF_dS = s.dMassFlux_dS[i];

      // Arithmetic average of fractional flow between local and neighbor (if any)
      real64 f_nei = f_loc;
      real64 df_nei_dP = 0.0;
      real64 df_nei_dS = 0.0;
      bool const hasNeighbor = (ner >= 0);
      if( hasNeighbor )
      {
        // compute neighbor fractional flow and its derivatives
        fracFlow( ner, nesr, nei, m_indep, f_nei, df_nei_dP, df_nei_dS );
      }

      // Upwind convex combination: beta = 1 if F >= 0 (use local), else 0 (use neighbor)
      real64 const beta = ( F >= 0.0 ) ? 1.0 : 0.0;
      real64 const f_int = beta * f_loc + ( 1.0 - beta ) * f_nei;

      // residual contribution and local Jacobians
      s.divSatFluxes += m_dt * F * f_int;
      // d/dP (local): dF/dP * f_int + F * beta * df_loc/dP (neighbor f has no local P dependence)
      s.dDivSatFluxes_dP += m_dt * ( dF_dP * f_int + F * beta * df_loc_dP );
      // d/dS (local): dF/dS * f_int + F * beta * df_loc/dS (neighbor f has no local S dependence)
      s.dDivSatFluxes_dS += m_dt * ( dF_dS * f_int + F * beta * df_loc_dS );
      // face pressure derivatives: only via F
      for( integer j=0; j<NUM_FACE; ++j ) s.dDivSatFluxes_dFaceVars[j] += m_dt * ( s.dMassFlux_dFacePres[i][j] * f_int );

      // neighbor Jacobian contributions on neighbor columns when applicable
      if( hasNeighbor )
      {
        // Neighbor global dof indices: pressure and saturation of independent phase
        globalIndex const neiP = m_elemDofNumber[ner][nesr][nei];
        globalIndex const neiS = neiP + 1;
        // Always append entries (values are be zero if beta == 1)
        s.neiCols[s.numNeiCols] = neiP;
        s.neiVals[s.numNeiCols] += m_dt * F * ( 1.0 - beta ) * df_nei_dP;
        ++s.numNeiCols;
        s.neiCols[s.numNeiCols] = neiS;
        s.neiVals[s.numNeiCols] += m_dt * F * ( 1.0 - beta ) * df_nei_dS;
        ++s.numNeiCols;
      }
    }
  }

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei, StackVariables & s, FUNC && ) const
  {
    if( m_elemGhostRank[ei] < 0 ){
      
      real64 const perm[3] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };
      
      // pressure equation
      IP::template compute< NUM_FACE >( m_nodePosition, m_transMultiplier, m_faceToNodes, m_elemToFaces[ei], m_elemCenter[ei], m_elemVolume[ei], perm, m_lengthTolerance, s.transMatrix );
      computeOverallMassFlux( ei, s );
      computeOverallMassFluxDivergence( ei, s );
      
      // saturation equation
      computeSaturationTransport( ei, s );
    }
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei, StackVariables & s ) const
  {
    if( m_elemGhostRank[ei] < 0 )
    {
      // Pressure equation row (cell pressure -> row 0 of cell block)
      localIndex const pressRow = s.cellRow;
      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[pressRow], s.divMassFluxes );
      real64 jacElemP[2] = { s.dDivMassFluxes_dP, s.dDivMassFluxes_dS };
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( pressRow, &s.elemCols[0], &jacElemP[0], 2 );
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( pressRow, &s.faceCols[0], &s.dDivMassFluxes_dFaceVars[0], NUM_FACE );
      
      // Saturation equation row (cell saturation -> row 1 of cell block)
      localIndex const satRow = s.cellRow + 1;
      RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[satRow], s.divSatFluxes );
      real64 jacElemS[2] = { s.dDivSatFluxes_dP, s.dDivSatFluxes_dS };
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( satRow, &s.elemCols[0], &jacElemS[0], 2 );
      m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( satRow, &s.faceCols[0], &s.dDivSatFluxes_dFaceVars[0], NUM_FACE );
      if( s.numNeiCols > 0 ){
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( satRow, &s.neiCols[0], &s.neiVals[0], s.numNeiCols );
      }
    }
    // face constraints unchanged: enforce hybrid face-pressure constraints using mass flux
    globalIndex const elemCol = s.elemCols[0];
    for( integer i=0; i<NUM_FACE; ++i )
    {
      if( m_faceGhostRank[m_elemToFaces[ei][i]] < 0 )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[s.faceRow[i]], s.MassFlux[i] );
        m_localMatrix.addToRow< parallelDeviceAtomic >( s.faceRow[i], &elemCol, &s.dMassFlux_dPres[i], 1 );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( s.faceRow[i], &s.faceCols[0], s.dMassFlux_dFacePres[i], NUM_FACE );
      }
    }
  }
private:
  globalIndex const m_rankOffset; localIndex const m_er; localIndex const m_esr; real64 const m_lengthTolerance; real64 const m_dt;
  arrayView1d< integer const > const m_elemGhostRank; ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumber;
  ElementRegionManager::ElementViewConst< arrayView2d< real64 const > > const m_phaseVolFracAll;
  arrayView1d< integer const > const m_faceGhostRank; arrayView1d< globalIndex const > const m_faceDofNumber;
  arrayView2d< localIndex const > const m_elemToFaces; arrayView2d< real64 const > const m_elemCenter; arrayView1d< real64 const > const m_elemVolume; arrayView1d< real64 const > const m_elemGravCoef;
  ArrayOfArraysView< localIndex const > const m_faceToNodes; arrayView1d< real64 const > const m_faceGravCoef; SortedArrayView< localIndex const > const m_regionFilter;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition; arrayView2d< localIndex const > const m_elemRegionList; arrayView2d< localIndex const > const m_elemSubRegionList; arrayView2d< localIndex const > const m_elemList;
  arrayView3d< real64 const > const m_elemPerm; arrayView1d< real64 const > const m_transMultiplier;
  arrayView1d< real64 const > const m_elemPres; arrayView1d< real64 const > const m_facePres;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens; arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;
  ElementRegionManager::ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseMobAll;
  ElementRegionManager::ElementViewConst< arrayView3d< real64 const, immiscibleFlow::USD_PHASE_DS > > const m_dPhaseMobAll;
  bool const m_assembleCellEquation;
  integer const m_indep;
  CRSMatrixView< real64, globalIndex const > const m_localMatrix; arrayView1d< real64 > const m_localRhs;
};

// Free function launcher for ElementBasedAssemblyKernel (avoids needing a static member)
// (Fixed: was using undefined GEOS_HOST_DEVICE_INLINE macro; use standard pattern instead.)
template< typename POLICY, integer NF, typename IPType >
GEOS_HOST_DEVICE inline void launchElementBasedAssemblyKernel( localIndex const n,
                                                               ElementBasedAssemblyKernel< NF, IPType > const k )
{
  forAll< POLICY >( n, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    typename ElementBasedAssemblyKernel< NF, IPType >::StackVariables s;
    k.setup( ei, s );
    k.compute( ei, s, NoOpFunc{} );
    k.complete( ei, s );
  } );
}

class ElementBasedAssemblyKernelFactory
{
public:
  template< typename POLICY >
  static void createAndLaunch( globalIndex const rankOffset,
                               localIndex const er,
                               localIndex const esr,
                               real64 const lengthTolerance,
                               string const elemDofKey,
                               string const faceDofKey,
                               string const solverName,
                               NodeManager const & nodeManager,
                               FaceManager const & faceManager,
                               ElementRegionManager const & elemManager,
                               CellElementSubRegion const & subRegion,
                               mimeticInnerProduct::MimeticInnerProductBase const & ipBase,
                               constitutive::TwoPhaseImmiscibleFluid const & fluid,
                               constitutive::PermeabilityBase const & permeability,
                               SortedArrayView< localIndex const > const & regionFilter,
                               integer const indepPhaseIndex,
                               real64 const & dt,
                               bool const assembleCellEq,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    mimeticInnerProductDispatch( ipBase, [&] ( auto const ip )
    {
      using IPType = TYPEOFREF( ip );
      internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
      {
        // persistent accessors
        auto dofNumberAccessor = elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
        dofNumberAccessor.setName( solverName + "/accessors/" + elemDofKey );
        auto satAccessor = elemManager.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );
        satAccessor.setName( solverName + "/accessors/" + fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );
        auto mobAccessor = elemManager.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseMobility::key() );
        mobAccessor.setName( solverName + "/accessors/" + fields::immiscibleMultiphaseFlow::phaseMobility::key() );
        auto dMobAccessor = elemManager.constructArrayViewAccessor< real64, 3 >( fields::immiscibleMultiphaseFlow::dPhaseMobility::key() );
        dMobAccessor.setName( solverName + "/accessors/" + fields::immiscibleMultiphaseFlow::dPhaseMobility::key() );
        ElementBasedAssemblyKernel< NF, IPType > k( rankOffset, er, esr, lengthTolerance,
                                                    faceDofKey, nodeManager, faceManager,
                                                    subRegion, dofNumberAccessor, satAccessor, fluid, permeability,
                                                    regionFilter, indepPhaseIndex, dt, assembleCellEq, localMatrix, localRhs,
                                                    mobAccessor, dMobAccessor );
        launchElementBasedAssemblyKernel< POLICY, NF, IPType >( subRegion.size(), k );
      } );
    } );
  }
};

/**
 * Accumulation kernel for MFD immiscible two-phase system.
 * Builds two equations per cell:
 *  - Row 0: total mass accumulation = PV * sum_i( rho_i * s_i ) - sum_i( phaseMass_n[i] )
 *  - Row 1: independent phase mass accumulation = PV * rho_indep * s_indep - phaseMass_n[indep]
 * with PV = volume * porosity. Jacobian wrt [pressure, s_indep].
 */
class AccumulationMFDKernel
{
public:
  static constexpr integer numEqn = 2; // total mass and independent phase mass
  static constexpr integer numDof = 2; // pressure and s_indep

  using Deriv = immiscibleFlow::DerivativeOffset;

  AccumulationMFDKernel( integer const indep,
                         globalIndex const rankOffset,
                         string const dofKey,
                         ElementSubRegionBase const & subRegion,
                         constitutive::TwoPhaseImmiscibleFluid const & fluid,
                         constitutive::CoupledSolidBase const & solid,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs )
    : m_indep( indep ),
      m_rankOffset( rankOffset ),
      m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
      m_elemGhostRank( subRegion.ghostRank() ),
      m_volume( subRegion.getElementVolume() ),
      m_porosity( solid.getPorosity() ),
      m_dPoro_dPres( solid.getDporosity_dPressure() ),
      m_phaseVolFrac( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >() ),
      m_phaseDens( fluid.phaseDensity() ),
      m_dPhaseDens( fluid.dPhaseDensity() ),
      m_phaseMass_n( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass_n >() ),
      m_localMatrix( localMatrix ),
      m_localRhs( localRhs )
  {}

  struct StackVariables
  {
    real64 PV = 0.0;
    real64 dPV_dP = 0.0;
    localIndex localRow = -1;
    globalIndex dofIndices[numDof]{};
    real64 localResidual[numEqn]{};
    real64 localJacobian[numEqn][numDof]{};
  };

  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const { return m_elemGhostRank[ei]; }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei, StackVariables & s ) const
  {
    s.PV = m_volume[ei] * m_porosity[ei][0];
    s.dPV_dP = m_volume[ei] * m_dPoro_dPres[ei][0];
    s.localRow = m_dofNumber[ei] - m_rankOffset;
    s.dofIndices[0] = m_dofNumber[ei];     // pressure
    s.dofIndices[1] = m_dofNumber[ei] + 1; // saturation of independent phase
  }

  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei, StackVariables & s ) const
  {
    // Current state
    real64 const s_ind = m_phaseVolFrac[ei][m_indep];
    integer const dep = 1 - m_indep;
    real64 const s_dep = m_phaseVolFrac[ei][dep]; // should already satisfy s_dep = 1 - s_ind

    real64 const rho_ind = m_phaseDens[ei][0][m_indep];
    real64 const rho_dep = m_phaseDens[ei][0][dep];
    real64 const drho_ind_dP = m_dPhaseDens[ei][0][m_indep][Deriv::dP];
    real64 const drho_dep_dP = m_dPhaseDens[ei][0][dep][Deriv::dP];

    // Previous masses
    real64 mass_n_ind = m_phaseMass_n[ei][m_indep];
    real64 mass_n_dep = m_phaseMass_n[ei][dep];
    // Fallback for first step: if both are zero, assume uninitialized and compute from current state
    if( mass_n_ind == 0.0 && mass_n_dep == 0.0 )
    {
      mass_n_ind = s.PV * rho_ind * s_ind;
      mass_n_dep = s.PV * rho_dep * s_dep;
    }

    // Equation 0: total mass
    real64 const mass_tot = s.PV * ( rho_ind * s_ind + rho_dep * s_dep );
    real64 const mass_tot_n = mass_n_ind + mass_n_dep;
    s.localResidual[0] += mass_tot - mass_tot_n;
    // d/dP
    real64 const sum_rho_s = rho_ind * s_ind + rho_dep * s_dep;
    real64 const sum_s_drho_dP = s_ind * drho_ind_dP + s_dep * drho_dep_dP;
    s.localJacobian[0][0] += s.dPV_dP * sum_rho_s + s.PV * sum_s_drho_dP;
    // d/dS_ind
    // d(sum_i rho_i s_i)/dS_ind = rho_ind - rho_dep (since s_dep = 1 - s_ind)
    s.localJacobian[0][1] += s.PV * ( rho_ind - rho_dep );

    // Equation 1: independent phase mass
    real64 const mass_ind = s.PV * rho_ind * s_ind;
    s.localResidual[1] += mass_ind - mass_n_ind;
    // d/dP
    s.localJacobian[1][0] += s.dPV_dP * ( rho_ind * s_ind ) + s.PV * ( s_ind * drho_ind_dP );
    // d/dS_ind
    s.localJacobian[1][1] += s.PV * rho_ind;
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const /*ei*/, StackVariables & s ) const
  {
    for( integer r = 0; r < numEqn; ++r )
    {
      m_localRhs[s.localRow + r] += s.localResidual[r];
      m_localMatrix.addToRow< serialAtomic >( s.localRow + r, s.dofIndices, s.localJacobian[r], numDof );
    }
  }

  template< typename POLICY, typename KERNEL_TYPE >
  static void launch( localIndex const n, KERNEL_TYPE const & k )
  {
    forAll< POLICY >( n, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( k.elemGhostRank( ei ) >= 0 ) return;
      typename KERNEL_TYPE::StackVariables s;
      k.setup( ei, s );
      k.computeAccumulation( ei, s );
      k.complete( ei, s );
    } );
  }

private:
  integer const m_indep;
  globalIndex const m_rankOffset;
  arrayView1d< globalIndex const > const m_dofNumber;
  arrayView1d< integer const > const m_elemGhostRank;
  arrayView1d< real64 const > const m_volume;
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseVolFrac;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseMass_n;
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  arrayView1d< real64 > const m_localRhs;
};

class AccumulationMFDKernelFactory
{
public:
  template< typename POLICY >
  static void createAndLaunch( globalIndex const rankOffset,
                               integer const indep,
                               string const dofKey,
                               ElementSubRegionBase const & subRegion,
                               constitutive::TwoPhaseImmiscibleFluid const & fluid,
                               constitutive::CoupledSolidBase const & solid,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    AccumulationMFDKernel kernel( indep, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    AccumulationMFDKernel::template launch< POLICY >( subRegion.size(), kernel );
  }
};

} // namespace immiscibleMultiphaseMFDKernels
} // namespace geos

#endif
