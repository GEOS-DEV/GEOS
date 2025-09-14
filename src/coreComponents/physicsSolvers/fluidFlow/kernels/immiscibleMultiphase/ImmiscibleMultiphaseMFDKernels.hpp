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
    
    // Keep pressure and saturation derivatives since lambda depends on both
    real64 Lambda = 0.0;   // == massLambda (before sign mapping)
    real64 dLambda_dP = 0.0;
    real64 dLambda_dS = 0.0; // derivative w.r.t. independent saturation (apply sign map below)
    
    // two-phase for now
    for( integer ip=0; ip<2; ++ip )
    {
      real64 const rho = m_phaseDens[ei][0][ip];
      real64 const drho_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
      real64 const lambda = m_phaseMobAll[m_er][m_esr][ei][ip];
      real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dP];
      real64 dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dS];

      // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
      real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
      real64 const dlambda_dS = sgnS * dlambda_dS_raw;
      
      // accumulate total mobility derivatives for rho_mix denominator
      Lambda += rho * lambda;
      dLambda_dP += drho_dP * lambda + rho * dlambda_dP;
      dLambda_dS += rho * dlambda_dS;
    }
    
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

        real64 const T_ij = m_dt * s.transMatrix[i][j];
        s.BuoyantFlux[i] += Lambda * T_ij * gravTerm;
        s.dBuoyantFlux_dPres[i] += T_ij * ( Lambda * dGravTerm_dP + dLambda_dP * gravTerm );
        s.dBuoyantFlux_dS[i] += T_ij * (dLambda_dS * gravTerm );
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
      real64 const rho = m_phaseDens[ei][0][ip];
      real64 const drho_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
      real64 const lambda = m_phaseMobAll[m_er][m_esr][ei][ip];
      real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dP];
      real64 dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei][ip][Deriv::dS];

      // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
      real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
      real64 const dlambda_dS = sgnS * dlambda_dS_raw;
      
      // accumulate total mobility derivatives for rho_mix denominator
      Lambda += rho * lambda;
      dLambda_dP += drho_dP * lambda + rho * dlambda_dP;
      dLambda_dS += rho * dlambda_dS;
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

      rho_hat += rho * (rho * lambda / Lambda);
      // numerator for d(rho_hat)/dP
      real64 const num_dP = Lambda * ( 2.0 * rho * lambda * drho_dP
                                + rho * rho * dlambda_dP )
                      - rho * rho * lambda * dLambda_dP;

      // numerator for d(rho_hat)/dS
      real64 const num_dS = Lambda * ( 2.0 * rho * lambda * drho_dS
                                + rho * rho * dlambda_dS )
                      - rho * rho * lambda * dLambda_dS;
      
      drho_hat_dP += num_dP / (Lambda * Lambda);
      drho_hat_dS += num_dS / (Lambda * Lambda);
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

        real64 const T_ij = m_dt * s.transMatrix[i][j];
        s.MassFlux[i] += Lambda * T_ij * potDif;
        s.dMassFlux_dPres[i] += T_ij * ( Lambda * dPotDif_dP + dLambda_dP * potDif );
        s.dMassFlux_dS[i] += T_ij * ( Lambda * dPotDif_dS + dLambda_dS * potDif );
        s.dMassFlux_dFacePres[i][j] += T_ij * ( Lambda * dPotDif_dFaceP );
      }
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
      s.divMassFluxes += F;
      // jacobians wrt element DOFs
      s.dDivMassFluxes_dP += dF_dP;
      s.dDivMassFluxes_dS += dF_dS;
      // wrt face pressures
      for( integer j=0; j<NUM_FACE; ++j ){
        s.dDivMassFluxes_dFaceVars[j] += s.dMassFlux_dFacePres[i][j];
      }
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
      integer const dep = 1 - indep;
      
      
      // Keep pressure and saturation derivatives since lambda depends on both
      real64 Lambda = 0.0;   // == massLambda (before sign mapping)
      real64 dLambda_dP = 0.0;
      real64 dLambda_dS = 0.0; // derivative w.r.t. independent saturation (apply sign map below)
      
      // two-phase for now
      for( integer ip=0; ip<2; ++ip )
      {
        real64 const rho = m_phaseDens[ei_local][0][ip];
        real64 const drho_dP = m_dPhaseDens[ei_local][0][ip][Deriv::dP];
        real64 const lambda = m_phaseMobAll[m_er][m_esr][ei_local][ip];
        real64 const dlambda_dP = m_dPhaseMobAll[m_er][m_esr][ei_local][ip][Deriv::dP];
        real64 dlambda_dS_raw = m_dPhaseMobAll[m_er][m_esr][ei_local][ip][Deriv::dS];

        // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
        real64 const sgnS = ( m_indep == 0 ? 1.0 : -1.0 );
        real64 const dlambda_dS = sgnS * dlambda_dS_raw;
        
        // accumulate total mobility derivatives for rho_mix denominator
        Lambda += rho * lambda;
        dLambda_dP += drho_dP * lambda + rho * dlambda_dP;
        dLambda_dS += rho * dlambda_dS;
      }
      
      // Phase mobilities
      real64 const lambda = m_phaseMobAll[er][esr][ei_local][indep];
      real64 const dlambda_dP = m_dPhaseMobAll[er][esr][ei_local][indep][Deriv::dP];
      // Map derivatives to configured independent saturation: if indep=1, d/dS1 = - d/dS0
      real64 const sgnS = ( indep == 0 ? 1.0 : -1.0 );
      real64 const dlambda_dS = sgnS * m_dPhaseMobAll[er][esr][ei_local][indep][Deriv::dS];

      // Phase densities (q=0 point) and pressure derivatives
      real64 const rho = m_phaseDens[ei_local][0][indep];
      real64 const drho_dP = m_dPhaseDens[ei_local][0][indep][Deriv::dP];
      real64 const drho_dS = 0;


      // df/dP
      real64 const num_dP = Lambda * (lambda * drho_dP + rho * dlambda_dP)
                    - rho * lambda * dLambda_dP;
      // df/dS
      real64 const num_dS = Lambda * (lambda * drho_dS + rho * dlambda_dS)
                    - rho * lambda * dLambda_dS;

      // Fractional flow (mass based)
      f = rho * lambda/ Lambda;
      df_dP = num_dP / (Lambda * Lambda);
      df_dS = num_dS / (Lambda * Lambda);
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
      localIndex ner = -1, nesr = -1, nei = -1;
      if( er0 == m_er && esr0 == m_esr && ei0 == ei ) { ner = er1; nesr = esr1; nei = ei1; }
      else if( er1 == m_er && esr1 == m_esr && ei1 == ei ) { ner = er0; nesr = esr0; nei = ei0; }

      // mass flux and derivatives for this one-sided face i (from cell to face)
      real64 const F = s.MassFlux[i];
      real64 const dF_dP = s.dMassFlux_dPres[i];
      real64 const dF_dS = s.dMassFlux_dS[i];
      
      // buoyant flux and derivatives for this one-sided face i (from cell to face)
      real64 const B = s.BuoyantFlux[i];
      real64 const dB_dP = s.dBuoyantFlux_dPres[i];
      real64 const dB_dS = s.dBuoyantFlux_dS[i];

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
      s.divSatFluxes += F * f_int;
      // d/dP (local): dF/dP * f_int + F * beta * df_loc/dP (neighbor f has no local P dependence)
      s.dDivSatFluxes_dP += ( dF_dP * f_int + F * beta * df_loc_dP );
      // d/dS (local): dF/dS * f_int + F * beta * df_loc/dS (neighbor f has no local S dependence)
      s.dDivSatFluxes_dS += ( dF_dS * f_int + F * beta * df_loc_dS );
      // face pressure derivatives
      for( integer j=0; j<NUM_FACE; ++j ){
        s.dDivSatFluxes_dFaceVars[j] += ( s.dMassFlux_dFacePres[i][j] * f_int );
      }

      // neighbor Jacobian contributions on neighbor columns when applicable
      if( hasNeighbor )
      {
        
        // Neighbor global dof indices: pressure and saturation of independent phase
        globalIndex const neiP = m_elemDofNumber[ner][nesr][nei];
        globalIndex const neiS = neiP + 1;
        
        // Always append entries (values are be zero if beta == 1)
        s.neiCols[s.numNeiCols] = neiP;
        s.neiVals[s.numNeiCols] += F * ( 1.0 - beta ) * df_nei_dP;
        ++s.numNeiCols;
        s.neiCols[s.numNeiCols] = neiS;
        s.neiVals[s.numNeiCols+1] += F * ( 1.0 - beta ) * df_nei_dS;
        ++s.numNeiCols;
        
        // perfom buoyancy contributions on internal facets only
        
        // Precompute local f and derivatives
        real64 f_dep_loc = 0.0, df_dep_loc_dP = 0.0, df_dep_loc_dS = 0.0;
        real64 f_dep_nei = 0.0, df_dep_nei_dP = 0.0, df_dep_nei_dS = 0.0;
        fracFlow( m_er, m_esr, ei, 1-m_indep, f_dep_loc, df_dep_loc_dP, df_dep_loc_dS );
        fracFlow( ner, nesr, nei, 1-m_indep, f_dep_nei, df_dep_nei_dP, df_dep_nei_dS );
        
        // Upwind convex combination: beta = 1 if F >= 0 (use local), else 0 (use neighbor)
        real64 const beta_g = ( B >= 0.0 ) ? 1.0 : 0.0;
        
        real64 const f_indep = beta_g * f_loc + ( 1.0 - beta_g ) * f_nei;
        real64 const df_indep_dP = beta_g * df_loc_dP + ( 1.0 - beta_g ) * df_nei_dP;
        real64 const df_indep_dS = beta_g * df_loc_dS + ( 1.0 - beta_g ) * df_nei_dS;
        
        real64 const f_dep = ( 1.0 - beta_g ) * f_dep_loc + beta_g * f_dep_nei;
        real64 const df_dep_dP = ( 1.0 - beta_g ) * df_dep_loc_dP + beta_g * df_dep_nei_dP;
        real64 const df_dep_dS = ( 1.0 - beta_g ) * df_dep_loc_dS + beta_g * df_dep_nei_dS;

        // residual contribution and local Jacobians
        s.divSatFluxes += f_indep * f_dep * B;
        // d/dP (local): df_indep_dP * f_dep * B + f_indep * df_dep_dP * B + f_indep * f_dep * dB_dP
        s.dDivSatFluxes_dP += ( df_indep_dP * f_dep * B + f_indep * df_dep_dP * B + f_indep * f_dep * dB_dP );
        // d/dS (local): df_indep_dS * f_dep * B + f_indep * df_dep_dS * B + f_indep * f_dep * dB_dS
        s.dDivSatFluxes_dS += ( df_indep_dS * f_dep * B + f_indep * df_dep_dS * B + f_indep * f_dep * dB_dS );
        
        s.neiVals[s.numNeiCols] += df_indep_dP * f_dep * B + f_indep * df_dep_dP * B + f_indep * f_dep * dB_dP;
        s.neiVals[s.numNeiCols+1] += df_indep_dS * f_dep * B + f_indep * df_dep_dS * B + f_indep * f_dep * dB_dS;
        
      }
    }
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
