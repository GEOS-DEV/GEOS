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
 * @file TemplatedFlowSolvers.cpp
 * @brief Implementation of template-based flow solver hierarchy
 */

#include "TemplatedFlowSolvers.hpp"
#include "FlowSolverBaseFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

// ========================================
// Base TemplatedFlowSolver Implementation
// ========================================

template< int NUM_PHASES, EquationType EQUATION_TYPE, bool HAS_CAPILLARITY, bool HAS_THERMAL >
TemplatedFlowSolver< NUM_PHASES, EQUATION_TYPE, HAS_CAPILLARITY, HAS_THERMAL >::
TemplatedFlowSolver( string const & name, Group * const parent )
  : FlowSolverBase( name, parent )
{
  // Set thermal flag if template parameter indicates thermal effects
  if( HAS_THERMAL )
  {
    m_isThermal = 1;
  }

  // Register phase-specific fields based on NUM_PHASES
  if constexpr ( NUM_PHASES > 1 )
  {
    // Register saturation fields for multiphase flow
    for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
    {
      string saturationKey = "saturation_phase" + std::to_string( phaseIndex );
      string saturationKey_n = saturationKey + "_n";
      
      // These would need to be registered in a registerDataOnMesh override
    }
  }

  // Register flux fields
  this->registerWrapper( "totalMassFlux", &m_totalMassFlux ).
    setDescription( "Total mass flux across all phases" );

  if constexpr ( HAS_CAPILLARITY )
  {
    this->registerWrapper( "totalCapillaryFlux", &m_totalCapillaryFlux ).
      setDescription( "Total capillary flux across all phases" );
  }
}

// ========================================
// PressureEquationSolver Implementation
// ========================================

template< int NUM_PHASES, bool HAS_CAPILLARITY, bool HAS_THERMAL >
PressureEquationSolver< NUM_PHASES, HAS_CAPILLARITY, HAS_THERMAL >::
PressureEquationSolver( string const & name, Group * const parent )
  : Base( name, parent )
{
  // Pressure equation specific initialization for mimetic discretization
}

template< int NUM_PHASES, bool HAS_CAPILLARITY, bool HAS_THERMAL >
void PressureEquationSolver< NUM_PHASES, HAS_CAPILLARITY, HAS_THERMAL >::
assembleSystem( real64 const time,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Get the mesh and mimetic discretization
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    // Update phase properties before assembly
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      updatePhaseProperties( subRegion );
    } );

    // Assemble pressure equation using mimetic discretization
    assembleMimeticSystem( time, dt, mesh, regionNames, dofManager, localMatrix, localRhs );
  } );
}

template< int NUM_PHASES, bool HAS_CAPILLARITY, bool HAS_THERMAL >
void PressureEquationSolver< NUM_PHASES, HAS_CAPILLARITY, HAS_THERMAL >::
computeFluxContributions( ElementSubRegionBase & subRegion, localIndex const stencilSize )
{
  // Get pressure field
  arrayView1d< real64 const > const pressure = subRegion.getField< flow::pressure >();

  if constexpr ( NUM_PHASES == 1 )
  {
    // Single-phase flux computation using mimetic discretization
    computeSinglePhaseMimeticFlux( subRegion, pressure, stencilSize );
  }
  else
  {
    // Multi-phase flux computation using mimetic discretization
    computeMultiPhaseMimeticFlux( subRegion, pressure, stencilSize );
  }
}

template< int NUM_PHASES, bool HAS_CAPILLARITY, bool HAS_THERMAL >
void PressureEquationSolver< NUM_PHASES, HAS_CAPILLARITY, HAS_THERMAL >::
updatePhaseProperties( ElementSubRegionBase & subRegion )
{
  // Update porosity and permeability using base class functionality
  if constexpr ( std::is_same_v< ElementSubRegionBase, CellElementSubRegion > )
  {
    this->updatePorosityAndPermeability( static_cast< CellElementSubRegion & >( subRegion ) );
  }
  else if constexpr ( std::is_same_v< ElementSubRegionBase, SurfaceElementSubRegion > )
  {
    this->updatePorosityAndPermeability( static_cast< SurfaceElementSubRegion & >( subRegion ) );
  }

  // Update phase-specific properties
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // Update fluid properties based on current pressure and temperature
  arrayView1d< real64 const > const pressure = subRegion.getField< flow::pressure >();
  arrayView1d< real64 const > const temperature = subRegion.getField< flow::temperature >();

  constitutive::ConstitutivePassThru< MultiFluidBase >::execute( fluid, [&]( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
    {
      fluidWrapper.updateState( k, pressure[k], temperature[k] );
    } );
  } );
}

// ========================================
// TransportEquationSolver Implementation
// ========================================

template< int NUM_PHASES, int TRANSPORT_INDEX, bool HAS_CAPILLARITY, bool HAS_THERMAL >
TransportEquationSolver< NUM_PHASES, TRANSPORT_INDEX, HAS_CAPILLARITY, HAS_THERMAL >::
TransportEquationSolver( string const & name, Group * const parent )
  : Base( name, parent )
{
  // Transport equation specific initialization for mimetic discretization
}

template< int NUM_PHASES, int TRANSPORT_INDEX, bool HAS_CAPILLARITY, bool HAS_THERMAL >
void TransportEquationSolver< NUM_PHASES, TRANSPORT_INDEX, HAS_CAPILLARITY, HAS_THERMAL >::
assembleSystem( real64 const time,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Assemble transport equation for phase TRANSPORT_INDEX using mimetic discretization
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      // Get saturation field for this transport phase
      string saturationKey = "saturation_phase" + std::to_string( TRANSPORT_INDEX );
      arrayView1d< real64 const > const saturation = subRegion.getField< real64_array >( saturationKey );
      arrayView1d< real64 const > const saturation_n = subRegion.getField< real64_array >( saturationKey + "_n" );

      // Compute transport residual using mimetic discretization
      arrayView1d< real64 > residual( subRegion.size() );
      computeTransportResidual( subRegion, saturation, saturation_n, residual );

      // Add to global system using mimetic discretization structure
      // Implementation depends on the specific mimetic DOF layout
    } );
  } );
}

// ========================================
// Specializations
// ========================================

// Single-phase specialization
template< bool HAS_THERMAL >
SinglePhaseFlowSolver< HAS_THERMAL >::
SinglePhaseFlowSolver( string const & name, Group * const parent )
  : Base( name, parent )
{
  // Single-phase specific initialization for mimetic discretization
  this->m_numDofPerCell = 1; // Only pressure DOF
}

template< bool HAS_THERMAL >
void SinglePhaseFlowSolver< HAS_THERMAL >::
assembleSystem( real64 const time,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  // Call base pressure equation assembly using mimetic discretization
  Base::assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  // Add single-phase specific terms if needed
  if constexpr ( HAS_THERMAL )
  {
    // Add thermal contributions using mimetic discretization
    assembleThermalContributions( time, dt, domain, dofManager, localMatrix, localRhs );
  }
}

// Two-phase specialization
template< bool HAS_CAPILLARITY, bool HAS_THERMAL >
TwoPhaseFlowSolver< HAS_CAPILLARITY, HAS_THERMAL >::
TwoPhaseFlowSolver( string const & name, Group * const parent )
  : Base( name, parent )
{
  // Two-phase specific initialization for mimetic discretization
  this->m_numDofPerCell = 2; // Pressure + 1 saturation DOF

  // Create transport solver for the first phase
  m_transportSolver = std::make_unique< TransportSolver >( name + "_transport", this );
}

template< bool HAS_CAPILLARITY, bool HAS_THERMAL >
void TwoPhaseFlowSolver< HAS_CAPILLARITY, HAS_THERMAL >::
assembleSystem( real64 const time,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  // Assemble pressure equation using mimetic discretization
  Base::assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  // Assemble transport equation using mimetic discretization
  m_transportSolver->assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  // Add coupling terms between pressure and transport equations (mimetic-specific)
  assembleMimeticCouplingTerms( time, dt, domain, dofManager, localMatrix, localRhs );
}

// Three-phase specialization
template< bool HAS_CAPILLARITY, bool HAS_THERMAL >
ThreePhaseFlowSolver< HAS_CAPILLARITY, HAS_THERMAL >::
ThreePhaseFlowSolver( string const & name, Group * const parent )
  : Base( name, parent )
{
  // Three-phase specific initialization for mimetic discretization
  this->m_numDofPerCell = 3; // Pressure + 2 saturation DOFs

  // Create transport solvers for the first two phases
  m_transportSolver1 = std::make_unique< TransportSolver1 >( name + "_transport1", this );
  m_transportSolver2 = std::make_unique< TransportSolver2 >( name + "_transport2", this );
}

template< bool HAS_CAPILLARITY, bool HAS_THERMAL >
void ThreePhaseFlowSolver< HAS_CAPILLARITY, HAS_THERMAL >::
assembleSystem( real64 const time,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  // Assemble pressure equation using mimetic discretization
  Base::assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  // Assemble transport equations using mimetic discretization
  m_transportSolver1->assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );
  m_transportSolver2->assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  // Add coupling terms between all equations (mimetic-specific)
  assembleMimeticCouplingTerms( time, dt, domain, dofManager, localMatrix, localRhs );
}

// Explicit template instantiations for common use cases
template class SinglePhaseFlowSolver< false >;
template class SinglePhaseFlowSolver< true >;

template class TwoPhaseFlowSolver< false, false >;
template class TwoPhaseFlowSolver< true, false >;
template class TwoPhaseFlowSolver< false, true >;
template class TwoPhaseFlowSolver< true, true >;

template class ThreePhaseFlowSolver< false, false >;
template class ThreePhaseFlowSolver< true, false >;
template class ThreePhaseFlowSolver< false, true >;
template class ThreePhaseFlowSolver< true, true >;

} // namespace geos
