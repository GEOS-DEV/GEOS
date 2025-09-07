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
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              localIndex const er,
                              localIndex const esr,
                              real64 const & lengthTolerance,
                              string const faceDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              CellElementSubRegion const & subRegion,
                              DofNumberAccessor const & elemDofNumberAccessor,
                              constitutive::TwoPhaseImmiscibleFluid const & fluid,
                              constitutive::PermeabilityBase const & permeability,
                              SortedArrayView< localIndex const > const & regionFilter,
                              real64 const & dt,
                              bool const assembleCellEq,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : m_rankOffset( rankOffset ), m_er( er ), m_esr( esr ), m_lengthTolerance( lengthTolerance ), m_dt( dt ),
      m_elemGhostRank( subRegion.ghostRank() ),
      m_elemDofNumber( elemDofNumberAccessor.toNestedViewConst() ),
      m_faceGhostRank( faceManager.ghostRank() ),
      m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
      m_elemToFaces( subRegion.faceList().toViewConst() ),
      m_elemCenter( subRegion.getElementCenter() ),
      m_elemVolume( subRegion.getElementVolume() ),
      m_elemGravCoef( subRegion.getField< fields::flow::gravityCoefficient >() ),
      m_faceToNodes( faceManager.nodeList().toViewConst() ),
      m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
      m_regionFilter( regionFilter ),
      m_nodePosition( nodeManager.referencePosition() ),
      m_elemRegionList( faceManager.elementRegionList() ),
      m_elemSubRegionList( faceManager.elementSubRegionList() ),
      m_elemList( faceManager.elementList() ),
      m_elemPerm( permeability.permeability() ),
      m_transMultiplier( faceManager.getField< fields::flow::transMultiplier >() ),
      m_elemPres( subRegion.getField< fields::flow::pressure >() ),
      m_facePres( faceManager.getField< fields::flow::facePressure >() ),
      m_phaseDens( fluid.phaseDensity() ),
      m_dPhaseDens( fluid.dPhaseDensity() ),
      m_phaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >() ),
      m_dPhaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >() ),
      m_assembleCellEquation( assembleCellEq ),
      m_localMatrix( localMatrix ), m_localRhs( localRhs )
  {}

  struct StackVariables
  {
    GEOS_HOST_DEVICE StackVariables(): transMatrix( NUM_FACE, NUM_FACE ) {}
    stackArray2d< real64, NUM_FACE * NUM_FACE > transMatrix;
    real64 oneSidedVolFlux[NUM_FACE]{}; // total volumetric flux (no mobility)
    real64 dOneSidedVolFlux_dPres[NUM_FACE]{}; // derivative wrt element pressure
    real64 dOneSidedVolFlux_dFacePres[NUM_FACE][NUM_FACE]{}; // derivative wrt face pressures
    real64 divMassFluxes = 0; // accumulation for cell eqn
    real64 dDivMassFluxes_dElemVars[NUM_FACE+1]{}; // [0] elem pressure, rest upwind neighbor pressures
    real64 dDivMassFluxes_dFaceVars[NUM_FACE]{};
    localIndex cellRow = 0; localIndex faceRow[NUM_FACE]{}; globalIndex elemCols[NUM_FACE+1]{}; globalIndex faceCols[NUM_FACE]{};
  };

  GEOS_HOST_DEVICE void setup( localIndex const ei, StackVariables & s ) const
  {
    s.cellRow = m_elemDofNumber[m_er][m_esr][ei] - m_rankOffset; s.elemCols[0] = m_elemDofNumber[m_er][m_esr][ei];
    for( integer f=0; f<NUM_FACE; ++f ) { localIndex lf = m_elemToFaces[ei][f]; s.faceRow[f] = m_faceDofNumber[lf]-m_rankOffset; s.faceCols[f]=m_faceDofNumber[lf]; }
  }

  GEOS_HOST_DEVICE void computeGradient( localIndex const ei, StackVariables & s ) const
  {
    for( integer i=0; i<NUM_FACE; ++i )
    {
      for( integer j=0; j<NUM_FACE; ++j )
      {
        real64 const ccPres = m_elemPres[ei];
        real64 const fPres = m_facePres[m_elemToFaces[ei][j]];
        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][j]];
        // Mixture density: rho_mix = (sum_i rho_i * lambda_i) / (sum_i lambda_i)
        // and derivative wrt pressure via quotient rule
        real64 sumLambda = 0.0;
        real64 dSumLambda_dP = 0.0;
        real64 sumRhoLambda = 0.0;
        real64 dSumRhoLambda_dP = 0.0;
        // two-phase for now
        for( integer ip=0; ip<2; ++ip )
        {
          real64 const rho = m_phaseDens[ei][0][ip];
          real64 const drho_dP = m_dPhaseDens[ei][0][ip][DerivMob::dP];
          real64 const lambda = m_phaseMob[ei][ip];
          real64 const dlambda_dP = m_dPhaseMob[ei][ip][DerivMob::dP];
          sumLambda += lambda;
          dSumLambda_dP += dlambda_dP;
          sumRhoLambda += rho * lambda;
          dSumRhoLambda_dP += drho_dP * lambda + rho * dlambda_dP;
        }
        // this can be eliminated as totall mobility is always > 0
        real64 const eps = 1e-30;
        real64 const denom = (fabs( sumLambda ) > eps) ? sumLambda : (sumLambda >= 0.0 ? eps : -eps);
        real64 const densMix = sumRhoLambda / denom;
        real64 const dDensMix_dP = ( dSumRhoLambda_dP * denom - sumRhoLambda * dSumLambda_dP ) / ( denom * denom );
        // potential difference terms
        real64 const presDif = ccPres - fPres;
        real64 const gravCoefDif = ccGravCoef - fGravCoef;
        real64 const gravTerm = densMix * gravCoefDif;
        real64 const dGravTerm_dP = dDensMix_dP * gravCoefDif;
        real64 const potDif = presDif - gravTerm;
        real64 const dPotDif_dP = 1.0 - dGravTerm_dP;
        real64 const dPotDif_dFaceP = -1.0;
        real64 const T_ij = s.transMatrix[i][j];
        s.oneSidedVolFlux[i] += T_ij * potDif;
        s.dOneSidedVolFlux_dPres[i] += T_ij * dPotDif_dP;
        s.dOneSidedVolFlux_dFacePres[i][j] += T_ij * dPotDif_dFaceP;
      }
    }
  }

  GEOS_HOST_DEVICE void computeFluxDivergence( localIndex const ei, StackVariables & s ) const
  {
    for( integer i=0; i<NUM_FACE; ++i )
    {
      localIndex local[3] = { m_er, m_esr, ei };
      localIndex neighbor[3] = { m_er, m_esr, ei };
      bool const hasNeigh = hybridFVMKernels::CellConnectivity::isNeighborFound( local, i, m_elemRegionList, m_elemSubRegionList, m_elemList, m_regionFilter, m_elemToFaces[ei], neighbor );
      // upwind choose element producing positive flux direction
      bool const useLocal = ( s.oneSidedVolFlux[i] >= 0 ) || !hasNeigh;
      localIndex const erU = useLocal ? local[0] : neighbor[0];
      localIndex const esrU = useLocal ? local[1] : neighbor[1];
      localIndex const eiU = useLocal ? local[2] : neighbor[2];
      // total mobility (sum phases) at upwind element
      real64 mobTot = 0.0; real64 dMobTot_dP = 0.0;
      for( integer ip=0; ip<2; ++ip ) // currently 2 phases
      {
        mobTot += m_phaseMob[eiU][ip];
        dMobTot_dP += m_dPhaseMob[eiU][ip][DerivMob::dP];
      }
      real64 const dt_mobTot = m_dt * mobTot;
      s.divMassFluxes += dt_mobTot * s.oneSidedVolFlux[i];
      s.dDivMassFluxes_dElemVars[0] += m_dt * mobTot * s.dOneSidedVolFlux_dPres[i];
      s.dDivMassFluxes_dElemVars[i+1] = m_dt * dMobTot_dP * s.oneSidedVolFlux[i];
      for( integer j=0; j<NUM_FACE; ++j )
      {
        s.dDivMassFluxes_dFaceVars[j] += dt_mobTot * s.dOneSidedVolFlux_dFacePres[i][j];
      }
      s.elemCols[i+1] = m_elemDofNumber[erU][esrU][eiU];
    }
  }

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE void compute( localIndex const ei, StackVariables & s, FUNC && ) const
  {
    real64 const perm[3] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };
    IP::template compute< NUM_FACE >( m_nodePosition, m_transMultiplier, m_faceToNodes, m_elemToFaces[ei], m_elemCenter[ei], m_elemVolume[ei], perm, m_lengthTolerance, s.transMatrix );
    computeGradient( ei, s );
    if( m_elemGhostRank[ei] < 0 ) computeFluxDivergence( ei, s );
  }

  GEOS_HOST_DEVICE void complete( localIndex const ei, StackVariables & s ) const
  {
    if( m_assembleCellEquation && m_elemGhostRank[ei] < 0 )
    {
      m_localRhs[s.cellRow] += s.divMassFluxes;
      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( s.cellRow, &s.elemCols[0], &s.dDivMassFluxes_dElemVars[0], NUM_FACE+1 );
      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( s.cellRow, &s.faceCols[0], &s.dDivMassFluxes_dFaceVars[0], NUM_FACE );
    }
    // face constraints unchanged
    globalIndex const elemCol = s.elemCols[0];
    for( integer i=0; i<NUM_FACE; ++i )
    {
      if( m_faceGhostRank[m_elemToFaces[ei][i]] < 0 )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[s.faceRow[i]], s.oneSidedVolFlux[i] );
        m_localMatrix.addToRow< parallelDeviceAtomic >( s.faceRow[i], &elemCol, &s.dOneSidedVolFlux_dPres[i], 1 );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( s.faceRow[i], &s.faceCols[0], s.dOneSidedVolFlux_dFacePres[i], NUM_FACE );
      }
    }
  }
private:
  globalIndex const m_rankOffset; localIndex const m_er; localIndex const m_esr; real64 const m_lengthTolerance; real64 const m_dt;
  arrayView1d< integer const > const m_elemGhostRank; ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumber;
  arrayView1d< integer const > const m_faceGhostRank; arrayView1d< globalIndex const > const m_faceDofNumber;
  arrayView2d< localIndex const > const m_elemToFaces; arrayView2d< real64 const > const m_elemCenter; arrayView1d< real64 const > const m_elemVolume; arrayView1d< real64 const > const m_elemGravCoef;
  ArrayOfArraysView< localIndex const > const m_faceToNodes; arrayView1d< real64 const > const m_faceGravCoef; SortedArrayView< localIndex const > const m_regionFilter;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition; arrayView2d< localIndex const > const m_elemRegionList; arrayView2d< localIndex const > const m_elemSubRegionList; arrayView2d< localIndex const > const m_elemList;
  arrayView3d< real64 const > const m_elemPerm; arrayView1d< real64 const > const m_transMultiplier;
  arrayView1d< real64 const > const m_elemPres; arrayView1d< real64 const > const m_facePres;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens; arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseMob; arrayView3d< real64 const, immiscibleFlow::USD_PHASE_DS > const m_dPhaseMob;
  bool const m_assembleCellEquation;
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
        // persistent accessor for element dof numbers
        auto dofNumberAccessor = elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
        dofNumberAccessor.setName( solverName + "/accessors/" + elemDofKey );
        ElementBasedAssemblyKernel< NF, IPType > k( rankOffset, er, esr, lengthTolerance,
                                                    faceDofKey, nodeManager, faceManager,
                                                    subRegion, dofNumberAccessor, fluid, permeability,
                                                    regionFilter, dt, assembleCellEq, localMatrix, localRhs );
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
    real64 const mass_n_ind = m_phaseMass_n[ei][m_indep];
    real64 const mass_n_dep = m_phaseMass_n[ei][dep];

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
