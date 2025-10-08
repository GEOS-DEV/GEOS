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
    // Buoyancy driven fluxes (one-sided per face i)
    real64 BuoyantFlux[NUM_FACE]{};
    real64 dBuoyantFlux_dPres[NUM_FACE]{};
    real64 dBuoyantFlux_dS[NUM_FACE]{};
    // PhaseMass (PPU) fluxes (one-sided per face i)
    real64 PhaseMassFlux[NUM_FACE][2]{};
    // Saturation equation: divergence of upwinded transported saturation with mass flux
    real64 divSatFluxes = 0.0;
    real64 dDivSatFluxes_dP = 0.0;
    real64 dDivSatFluxes_dS = 0.0;
    real64 dDivSatFluxes_dFaceVars[NUM_FACE]{};
    // Neighbor pressure and saturation couplings (via upwind when inflow)
    localIndex numNeiCols = 0;
    globalIndex neiCols[2*NUM_FACE]{}; // allow both P and S per neighbor face
    real64 neiVals[2*NUM_FACE]{};
    // Row/col bookkeeping
    localIndex cellRow = 0;
    localIndex faceRow[NUM_FACE]{};
    globalIndex elemCols[2]{}; // [P, S_indep]
    globalIndex faceCols[NUM_FACE]{};
  };

  GEOS_HOST_DEVICE
  void setup( localIndex const ei, StackVariables & s ) const
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
  
  GEOS_HOST_DEVICE
  void computeBuoyancyDrivenFlux( localIndex const ei, StackVariables & s ) const
  {
    using Deriv = DerivMob; // alias for readability
    
    // Compute delta_rho and its pressure derivative
    integer const dep = 1 - m_indep;
    real64 const rho_dep = m_phaseDens[ei][0][dep];
    real64 const rho_indep = m_phaseDens[ei][0][m_indep];
    real64 const drho_dep_dP = m_dPhaseDens[ei][0][dep][Deriv::dP];
    real64 const drho_indep_dP = m_dPhaseDens[ei][0][m_indep][Deriv::dP];

    real64 const delta_rho = rho_dep - rho_indep;
    real64 const ddelta_rho_dP = drho_dep_dP - drho_indep_dP;
    
    for( integer i=0; i<NUM_FACE; ++i )
    {
      for( integer j=0; j<NUM_FACE; ++j )
      {
        // Local gravity terms (cell-centered and face values associated to j)
        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][j]];

        // Potential difference and its derivatives
        real64 const gravCoefDif = ccGravCoef - fGravCoef;
        real64 const gravTerm = delta_rho * gravCoefDif;
        real64 const dGravTerm_dP = ddelta_rho_dP * gravCoefDif;
        real64 const dGravTerm_dS = 0.0;

        real64 const T_ij = s.transMatrix[i][j];
        s.BuoyantFlux[i] += T_ij * gravTerm;
        s.dBuoyantFlux_dPres[i] += T_ij * (dGravTerm_dP);
        s.dBuoyantFlux_dS[i] += T_ij * (dGravTerm_dS);
      }
    }
  }
  
  GEOS_HOST_DEVICE
  void computeOverallMassFlux( localIndex const ei, StackVariables & s ) const
  {
    using Deriv = DerivMob; // alias for readability
    
    // Keep pressure and saturation derivatives since lambda depends on both
    real64 Lambda = 0.0;   // == massLambda (before sign mapping)
    real64 dLambda_dP = 0.0;
    real64 dLambda_dS = 0.0; // derivative w.r.t. independent saturation (apply sign map below)
    
    // two-phase for now
    for( integer ip=0; ip<2; ++ip )
    {
      // this mobility are mass mobilities (rho * kr / mu).
      real64 const lambda = m_phaseMobAll[m_er][m_esr][ei][ip];
      real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dP];
      real64 dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dS];

      // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
      real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
      real64 const dlambda_dS = sgnS * dlambda_dS_raw;
      
      // accumulate total mobility derivatives for rho_mix denominator
      Lambda += lambda;
      dLambda_dP += dlambda_dP;
      dLambda_dS += dlambda_dS;
    }
    
    // Mixture density rho_hat = (sum_i rho_i * (lambda_i / (sum_i lambda_i))
    real64 rho_hat = 0.0;
    real64 drho_hat_dP = 0.0;
    real64 drho_hat_dS = 0.0;
    for( integer ip=0; ip<2; ++ip )
    {
      real64 const rho = m_phaseDens[ei][0][ip];
      real64 const drho_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
      real64 const drho_dS = 0.0;
      real64 const lambda = m_phaseMobAll[m_er][m_esr][ei][ip];
      real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dP];
      real64 dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dS];

      // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
      real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
      real64 const dlambda_dS = sgnS * dlambda_dS_raw;

      rho_hat += rho * (lambda / Lambda);
      drho_hat_dP += (drho_dP * lambda * Lambda +
                      rho * (dlambda_dP * Lambda - lambda * dLambda_dP)) / (Lambda * Lambda);
      drho_hat_dS += (drho_dS * lambda * Lambda +
                      rho * (dlambda_dS * Lambda - lambda * dLambda_dS)) / (Lambda * Lambda);
    }
    
    for( integer i=0; i<NUM_FACE; ++i )
    {
      for( integer j=0; j<NUM_FACE; ++j )
      {
        // Local pressure and gravity terms (cell-centered and face values associated to j)
        real64 const ccPres = m_elemPres[ei];
        real64 const fPres = m_facePres[m_elemToFaces[ei][j]];
        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][j]];

        // Potential difference and its derivatives
        real64 const presDif = ccPres - fPres;
        real64 const gravCoefDif = ccGravCoef - fGravCoef;
        real64 const gravTerm = rho_hat * gravCoefDif;
        real64 const dGravTerm_dP = drho_hat_dP * gravCoefDif;
        real64 const dGravTerm_dS = drho_hat_dS * gravCoefDif;
        real64 const potDif = presDif - gravTerm;
        real64 const dPotDif_dP = 1.0 - dGravTerm_dP;
        real64 const dPotDif_dS = - dGravTerm_dS;
        real64 const dPotDif_dFaceP = -1.0;

        real64 const T_ij = s.transMatrix[i][j];
        s.MassFlux[i] += Lambda * T_ij * potDif;
        s.dMassFlux_dPres[i] += T_ij * ( Lambda * dPotDif_dP + dLambda_dP * potDif );
        s.dMassFlux_dS[i] += T_ij * ( Lambda * dPotDif_dS + dLambda_dS * potDif );
        s.dMassFlux_dFacePres[i][j] += T_ij * ( Lambda * dPotDif_dFaceP );
      }
    }
  }

  GEOS_HOST_DEVICE
  void computePhaseMassFlux( localIndex const ei, StackVariables & s ) const
  {
    // PhaseMass (PPU) fluxes for each phase (one-sided per face i)
    // References:
    // [1] Brenier, Y., & Jaffré, J. (1991). SIAM Journal on Numerical Analysis. DOI:10.1137/0728036
    // [2] Kwok, F., & Tchelepi, H. A. (2008). SIAM Journal on Numerical Analysis. DOI:10.1137/070703922
    integer const dep = 1 - m_indep;
    // Local (cell-centered) mass mobilities and fractional flows
    real64 const lambda_ind = m_phaseMobAll[m_er][m_esr][ei][m_indep];
    real64 const lambda_dep = m_phaseMobAll[m_er][m_esr][ei][dep];
    real64 const lambda_tot = lambda_ind + lambda_dep + 1e-32; // tiny guard
    real64 const f_ind = lambda_ind / lambda_tot;
    real64 const f_dep = lambda_dep / lambda_tot;
    real64 const l_sum = lambda_tot; // shorthand
    real64 const buoyancy_factor = f_ind * f_dep * l_sum; // matches earlier saturation treatment

    for( integer i=0; i<NUM_FACE; ++i )
    {
      real64 const F = s.MassFlux[i];      // total mass flux (one-sided cell->face)
      real64 const B = s.BuoyantFlux[i];   // buoyancy transmissibility (delta_rho * T * gravCoefDiff)
      // Phase indep
      s.PhaseMassFlux[i][m_indep] = f_ind * F + buoyancy_factor * B;
      // Phase dep
      s.PhaseMassFlux[i][dep]     = f_dep * F - buoyancy_factor * B;
    }
  }
  
  GEOS_HOST_DEVICE
  void computeOverallMassFluxDivergence( localIndex const ei, StackVariables & s ) const
  {
    GEOS_UNUSED_VAR( ei );
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
      for( integer j=0; j<NUM_FACE; ++j ){
        s.dDivMassFluxes_dFaceVars[j] += m_dt * s.dMassFlux_dFacePres[i][j];
      }
    }
  }

  GEOS_HOST_DEVICE
  void computeSaturationTransport(localIndex const ei, StackVariables & s) const
  {
    // Upwind options
    enum UpwindScheme : int { HU = 0, PPU = 1 };
    UpwindScheme const scheme = UpwindScheme::PPU; // expose later if needed

    using Deriv = DerivMob; // keep your existing derivative enum

    // --- Small types to hold eta and upwinded data (using loc/nei naming) ---
    struct EtaData {
      real64 valLoc;   real64 dP_Loc;   real64 dS_Loc;
      real64 valNei;   real64 dP_Nei;   real64 dS_Nei;
    };

    struct UpwindedEta {
      real64 alpha;    // upwind coefficient (1 -> loc, 0 -> nei)
      real64 value;
      real64 dP_Loc;   real64 dS_Loc;
      real64 dP_Nei;   real64 dS_Nei;
    };

    // ------------------------------
    // 1) mobility builder: eta = lambda_gamma for a given phase
    // ------------------------------
    auto buildMobility = [this] GEOS_HOST_DEVICE (localIndex er, localIndex esr, localIndex ei_local, integer phase) -> EtaData {
      EtaData out{};
      out.valLoc = m_phaseMobAll[er][esr][ei_local][phase];
      out.dP_Loc = m_dPhaseMobAll[er][esr][ei_local][phase][Deriv::dP];
      real64 const dlambda_dS_raw = m_dPhaseMobAll[er][esr][ei_local][phase][Deriv::dS];
      real64 const sgnS = (m_indep == 0 ? 1.0 : -1.0);
      out.dS_Loc = sgnS * dlambda_dS_raw;
      // initialize neighbor to equal local; caller will overwrite valNei/dNei if neighbor exists
      out.valNei = out.valLoc;
      out.dP_Nei = out.dP_Loc;
      out.dS_Nei = out.dS_Loc;
      return out;
    };

    // ------------------------------
    // 2) fractional flow builder: eta = f_gamma = lambda_gamma / Lambda
    //    uses buildMobility for each phase and computes derivatives for local cell
    // ------------------------------
    auto buildFracFlow = [buildMobility] GEOS_HOST_DEVICE (localIndex er, localIndex esr, localIndex ei_local, integer phase) -> EtaData {
      EtaData mob0 = buildMobility(er, esr, ei_local, 0);
      EtaData mob1 = buildMobility(er, esr, ei_local, 1);

      real64 const Lambda     = mob0.valLoc + mob1.valLoc;
      real64 const dLambda_dP = mob0.dP_Loc + mob1.dP_Loc;
      real64 const dLambda_dS = mob0.dS_Loc + mob1.dS_Loc;

      EtaData out{};
      EtaData const & target = (phase == 0 ? mob0 : mob1);

      real64 const lambda = target.valLoc;
      real64 const dl_dP  = target.dP_Loc;
      real64 const dl_dS  = target.dS_Loc;

      // value and local derivatives (neigh init to local)
      out.valLoc = lambda / Lambda;
      out.dP_Loc = ( dl_dP * Lambda - lambda * dLambda_dP ) / ( Lambda * Lambda );
      out.dS_Loc = ( dl_dS * Lambda - lambda * dLambda_dS ) / ( Lambda * Lambda );

      out.valNei = out.valLoc;
      out.dP_Nei = out.dP_Loc;
      out.dS_Nei = out.dS_Loc;
      return out;
    };

    // ------------------------------
    // 3) generic applyUpwind helper:
    //    given eta data and scalar beta, produce upwinded value and derivative partitions
    //    alpha = 1 if beta >= 0 (use loc), alpha = 0 otherwise (use nei)
    // ------------------------------
    auto applyUpwind = [] GEOS_HOST_DEVICE (EtaData const & eta, real64 beta) -> UpwindedEta {
      UpwindedEta u{};
      u.alpha = (beta >= 0.0 ? 1.0 : 0.0); // simple step upwinding (HU/PPU use different betas)
      u.value = u.alpha * eta.valLoc + (1.0 - u.alpha) * eta.valNei;
      u.dP_Loc = u.alpha * eta.dP_Loc;
      u.dS_Loc = u.alpha * eta.dS_Loc;
      u.dP_Nei = (1.0 - u.alpha) * eta.dP_Nei;
      u.dS_Nei = (1.0 - u.alpha) * eta.dS_Nei;
      return u;
    };

    // ------------------------------
    // 4) main face loop (almost identical structure but clearer naming)
    // ------------------------------
    for (integer i = 0; i < NUM_FACE; ++i)
    {
      // face / neighbor bookkeeping (loc cell is `ei`)
      localIndex const lf = m_elemToFaces[ei][i];
      localIndex const er0 = m_elemRegionList[lf][0];
      localIndex const esr0 = m_elemSubRegionList[lf][0];
      localIndex const ei0 = m_elemList[lf][0];
      localIndex const er1 = m_elemRegionList[lf][1];
      localIndex const esr1 = m_elemSubRegionList[lf][1];
      localIndex const ei1 = m_elemList[lf][1];

      localIndex ner = -1, nesr = -1, nei = -1;
      if (er0 == m_er && esr0 == m_esr && ei0 == ei) { ner = er1; nesr = esr1; nei = ei1; }
      else if (er1 == m_er && esr1 == m_esr && ei1 == ei) { ner = er0; nesr = esr0; nei = ei0; }
      bool const hasNeighbor = (ner >= 0);

      // fluxes and their derivatives
      real64 const F = s.MassFlux[i];
      real64 const dF_dP = s.dMassFlux_dPres[i];
      real64 const dF_dS = s.dMassFlux_dS[i];
      real64 const B = s.BuoyantFlux[i];
      real64 const dB_dP = s.dBuoyantFlux_dPres[i];
      real64 const dB_dS = s.dBuoyantFlux_dS[i];
      real64 const F_indep = s.PhaseMassFlux[i][m_indep];
      real64 const F_dep   = s.PhaseMassFlux[i][1 - m_indep];

      // build eta data for independent and dependent phases (loc values)
      EtaData f_indep = buildFracFlow(m_er, m_esr, ei, m_indep);
      EtaData f_dep   = buildFracFlow(m_er, m_esr, ei, 1 - m_indep);
      EtaData l_indep = buildMobility(m_er, m_esr, ei, m_indep);
      EtaData l_dep   = buildMobility(m_er, m_esr, ei, 1 - m_indep);

      // if neighbor exists, fetch neighbor loc-values and store them into .valNei/.d*_Nei
      if (hasNeighbor) {
        EtaData f_indep_nei = buildFracFlow(ner, nesr, nei, m_indep);
        f_indep.valNei = f_indep_nei.valLoc; f_indep.dP_Nei = f_indep_nei.dP_Loc; f_indep.dS_Nei = f_indep_nei.dS_Loc;

        EtaData f_dep_nei   = buildFracFlow(ner, nesr, nei, 1 - m_indep);
        f_dep.valNei = f_dep_nei.valLoc; f_dep.dP_Nei = f_dep_nei.dP_Loc; f_dep.dS_Nei = f_dep_nei.dS_Loc;

        EtaData l_indep_nei = buildMobility(ner, nesr, nei, m_indep);
        l_indep.valNei = l_indep_nei.valLoc; l_indep.dP_Nei = l_indep_nei.dP_Loc; l_indep.dS_Nei = l_indep_nei.dS_Loc;

        EtaData l_dep_nei   = buildMobility(ner, nesr, nei, 1 - m_indep);
        l_dep.valNei = l_dep_nei.valLoc; l_dep.dP_Nei = l_dep_nei.dP_Loc; l_dep.dS_Nei = l_dep_nei.dS_Loc;
      }

      // ----- Convective term (use chosen advecting scalar for convective upwinding) -----
      real64 const beta_conv = (scheme == UpwindScheme::PPU ? F_indep : F);
      UpwindedEta f_conv = applyUpwind(f_indep, beta_conv);

      s.divSatFluxes += m_dt * F * f_conv.value;
      s.dDivSatFluxes_dP += m_dt * (dF_dP * f_conv.value + F * f_conv.dP_Loc);
      s.dDivSatFluxes_dS += m_dt * (dF_dS * f_conv.value + F * f_conv.dS_Loc);

      // chain rule contributions to face-pres DOFs
      for (integer j = 0; j < NUM_FACE; ++j)
        s.dDivSatFluxes_dFaceVars[j] += m_dt * s.dMassFlux_dFacePres[i][j] * f_conv.value;

      // ----- Buoyancy terms (use different betas to choose upwind direction) -----
      real64 const beta_indep = (scheme == UpwindScheme::PPU ? F_indep : B);
      real64 const beta_dep   = (scheme == UpwindScheme::PPU ? F_dep   : -B);

      UpwindedEta f_up_indep = applyUpwind(f_indep, beta_indep);
      UpwindedEta f_up_dep   = applyUpwind(f_dep,   beta_dep);
      UpwindedEta l_up_indep = applyUpwind(l_indep, beta_indep);
      UpwindedEta l_up_dep   = applyUpwind(l_dep,   beta_dep);

      real64 const f_ind = f_up_indep.value;
      real64 const f_dep_v = f_up_dep.value;
      real64 const l_sum = l_up_indep.value + l_up_dep.value;
      real64 const buoyancy_factor = f_ind * f_dep_v * l_sum;

      s.divSatFluxes += m_dt * buoyancy_factor * B;

      // local contribution derivative helper (derivative wrt local P or S)
      auto dBuoy_dLocal = [&](bool dP) -> real64 {
        real64 df_ind = dP ? f_up_indep.dP_Loc : f_up_indep.dS_Loc;
        real64 df_dep = dP ? f_up_dep.dP_Loc   : f_up_dep.dS_Loc;
        real64 dl_ind = dP ? l_up_indep.dP_Loc  : l_up_indep.dS_Loc;
        real64 dl_dep = dP ? l_up_dep.dP_Loc    : l_up_dep.dS_Loc;
        real64 dB = dP ? dB_dP : dB_dS;
        // product rule: B * d(buoy_factor) + buoy_factor * dB
        return m_dt * ( B * ( (df_ind * f_dep_v + f_ind * df_dep) * l_sum + f_ind * f_dep_v * (dl_ind + dl_dep) ) + buoyancy_factor * dB );
      };

      s.dDivSatFluxes_dP += dBuoy_dLocal(true);
      s.dDivSatFluxes_dS += dBuoy_dLocal(false);

      // neighbor contributions (if neighbor exists) -- place in neighbor columns
      if (hasNeighbor) {
        // convective neighbor contribution
        real64 const conv_dP_nei = m_dt * F * f_conv.dP_Nei;
        real64 const conv_dS_nei = m_dt * F * f_conv.dS_Nei;

        // buoyancy neighbor contribution: derivative of buoy_factor with respect to neighbor P/S
        auto dBuoy_dNei = [&](bool dP) -> real64 {
          real64 df_ind = dP ? f_up_indep.dP_Nei : f_up_indep.dS_Nei;
          real64 df_dep = dP ? f_up_dep.dP_Nei   : f_up_dep.dS_Nei;
          real64 dl_ind = dP ? l_up_indep.dP_Nei  : l_up_indep.dS_Nei;
          real64 dl_dep = dP ? l_up_dep.dP_Nei    : l_up_dep.dS_Nei;
          return m_dt * B * ( (df_ind * f_dep_v + f_ind * df_dep) * l_sum + f_ind * f_dep_v * (dl_ind + dl_dep) );
        };

        real64 const buoy_dP_nei = dBuoy_dNei(true);
        real64 const buoy_dS_nei = dBuoy_dNei(false);

        globalIndex const neiP = m_elemDofNumber[ner][nesr][nei];
        globalIndex const neiS = neiP + 1;
        s.neiCols[s.numNeiCols] = neiP;
        s.neiVals[s.numNeiCols] += conv_dP_nei + buoy_dP_nei;
        s.numNeiCols += 1;
        s.neiCols[s.numNeiCols] = neiS;
        s.neiVals[s.numNeiCols] += conv_dS_nei + buoy_dS_nei;
        s.numNeiCols += 1;
      }
    } // end face loop
  }


  GEOS_HOST_DEVICE
  void compute( localIndex const ei, StackVariables & s) const
  {
    if( m_elemGhostRank[ei] < 0 ){
      
      real64 const perm[3] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };
      
      // pressure equation
      IP::template compute< NUM_FACE >( m_nodePosition, m_transMultiplier, m_faceToNodes, m_elemToFaces[ei], m_elemCenter[ei], m_elemVolume[ei], perm, m_lengthTolerance, s.transMatrix );
      computeOverallMassFlux( ei, s );
      computeOverallMassFluxDivergence( ei, s );
      
      // saturation equation
      computeBuoyancyDrivenFlux(ei, s);
      computePhaseMassFlux( ei, s );
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
    globalIndex const elemCol_p = s.elemCols[0];
    globalIndex const elemCol_s = s.elemCols[0] + 1;
    for( integer i=0; i<NUM_FACE; ++i )
    {
      if( m_faceGhostRank[m_elemToFaces[ei][i]] < 0 )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[s.faceRow[i]], s.MassFlux[i] );
        m_localMatrix.addToRow< parallelDeviceAtomic >( s.faceRow[i], &elemCol_p, &s.dMassFlux_dPres[i], 1 );
        m_localMatrix.addToRow< parallelDeviceAtomic >( s.faceRow[i], &elemCol_s, &s.dMassFlux_dS[i], 1 );
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
    k.compute( ei, s );
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
        ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor = elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
        dofNumberAccessor.setName( solverName + "/accessors/" + elemDofKey );
        ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > satAccessor = elemManager.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );
        satAccessor.setName( solverName + "/accessors/" + fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );
        ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > mobAccessor = elemManager.constructArrayViewAccessor< real64, 2 >( fields::immiscibleMultiphaseFlow::phaseMobility::key() );
        mobAccessor.setName( solverName + "/accessors/" + fields::immiscibleMultiphaseFlow::phaseMobility::key() );
        ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dMobAccessor = elemManager.constructArrayViewAccessor< real64, 3 >( fields::immiscibleMultiphaseFlow::dPhaseMobility::key() );
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
