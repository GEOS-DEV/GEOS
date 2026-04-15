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

#include "FieldApplicator.hpp"
#include "events/tasks/TasksManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationImpl.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/SurfaceElementSubRegion.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"

namespace geos
{
using namespace dataRepository;

FieldApplicator::FieldApplicator( const string & name,
                                  Group * const parent ):
  TaskBase( name, parent ),
  m_fieldSpecificationNames(),
  m_solverName(),
  m_targetRegions()
{

  registerWrapper( viewKeyStruct::fieldSpecificationNamesString(), &m_fieldSpecificationNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Array containing the field specifications to apply to newly created elements" );

  registerWrapper( viewKeyStruct::solverNameString(), &m_solverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the flow solver to use for fluid state initialization (required for compositional flow)" );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegions ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Target element regions to apply field specifications to. If empty, applies to all regions." );

}

FieldApplicator::~FieldApplicator() = default;

void FieldApplicator::postInputInitialization()
{
  TaskBase::postInputInitialization();
}


bool
FieldApplicator::
  execute( real64 const time_n,
           real64 const dt,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n, dt, cycleNumber, eventCounter, eventProgress );

  FieldSpecificationManager & fsm = FieldSpecificationManager::getInstance();

  // Check whether any of the field specifications is an EquilibriumInitialCondition.
  // If so, we delegate the full re-initialization to the flow solver, which correctly
  // handles multi-phase hydrostatic pressure computation and all derived state
  // (phaseVolumeFraction, phaseMobility, RelPerm_phaseRelPerm, fracPerm_permeability,
  //  fracPorosity_porosity, etc.)
  bool hasEquilibriumIC = false;
  for( string const & fsName : m_fieldSpecificationNames )
  {
    if( fsm.hasGroup( fsName ) )
    {
      FieldSpecification const & fs = fsm.getGroup< FieldSpecification >( fsName );
      if( dynamic_cast< EquilibriumInitialCondition const * >( &fs ) != nullptr )
      {
        hasEquilibriumIC = true;
        break;
      }
    }
  }

  if( hasEquilibriumIC )
  {
    // Find the flow solver to delegate initialization to.
    // Use m_solverName if provided, otherwise find the first FlowSolverBase.
    FlowSolverBase * flowSolver = nullptr;
    Group & solversGroup = this->getGroupByPath< Group >( "/Problem/Solvers" );

    if( !m_solverName.empty() )
    {
      flowSolver = solversGroup.getGroupPointer< FlowSolverBase >( m_solverName );
      GEOS_THROW_IF( flowSolver == nullptr,
                     GEOS_FMT( "FieldApplicator '{}': solver '{}' not found or is not a FlowSolverBase.",
                               getName(), m_solverName ),
                     InputError, getDataContext() );
    }
    else
    {
      solversGroup.forSubGroups< FlowSolverBase >( [&]( FlowSolverBase & solver )
      {
        if( flowSolver == nullptr )
        {
          flowSolver = &solver;
        }
      } );
    }

    GEOS_THROW_IF( flowSolver == nullptr,
                   GEOS_FMT( "FieldApplicator '{}': No FlowSolverBase found to delegate HydrostaticEquilibrium "
                             "initialization. Please specify solverName in the FieldApplicator XML tag.",
                             getName() ),
                   InputError, getDataContext() );

    GEOS_LOG_RANK_0( GEOS_FMT( "FieldApplicator '{}': delegating HydrostaticEquilibrium re-initialization "
                               "to solver '{}'.", getName(), flowSolver->getName() ) );

    // Re-run the solver's full state initialization, which internally calls (in order):
    //   1. computeHydrostaticEquilibrium  → correctly assigns pressure, temperature, and
    //                                       globalCompFraction to ALL subregions including the
    //                                       newly created fracture elements, using the full
    //                                       multi-phase hydrostatic pressure algorithm
    //   2. initializePorosityAndPermeability  → fracPorosity_porosity, fracPerm_permeability
    //   3. initializeHydraulicAperture        → hydraulicAperture for fractures
    //   4. initializeFluidState               → globalCompDensity, phaseVolumeFraction,
    //                                           RelPerm_phaseRelPerm, phaseMobility,
    //                                           fluid.saveConvergedState() (_n arrays)
    flowSolver->initializeState( domain );

    return false;
  }

  // For non-equilibrium field specifications, apply them element by element as before.
  for( string const & fsName : m_fieldSpecificationNames )
  {
    GEOS_THROW_IF( !fsm.hasGroup( fsName ),
                   GEOS_FMT( "{}: FieldSpecification named {} not found",
                             getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ),
                             fsName ),
                   InputError, getWrapperDataContext( viewKeyStruct::fieldSpecificationNamesString() ) );

    FieldSpecification const & fs = fsm.getGroup< FieldSpecification >( fsName );

    for( auto & meshBodyPair : domain.getMeshBodies().getSubGroups() )
    {
      if( MeshBody * const meshBody = dynamic_cast< MeshBody * >( meshBodyPair.second ) )
      {
        for( auto & meshLevelPair : meshBody->getMeshLevels().getSubGroups() )
        {
          if( MeshLevel * const meshLevel = dynamic_cast< MeshLevel * >( meshLevelPair.second ) )
          {
            FieldSpecificationImpl::apply< ElementSubRegionBase >( fs,
                                                                   *meshLevel,
                                                                   [&]( FieldSpecification const & bc,
                                                                        string const &,
                                                                        SortedArrayView< localIndex const > const & targetSet,
                                                                        ElementSubRegionBase & subRegion,
                                                                        string const fieldName )
            {
              // If targetRegions is specified, only apply to matching regions
              if( !m_targetRegions.empty() )
              {
                string const & regionName = subRegion.getParent().getName();
                bool isTargetRegion = false;
                for( string const & targetRegion : m_targetRegions )
                {
                  if( regionName == targetRegion )
                  {
                    isTargetRegion = true;
                    break;
                  }
                }
                if( !isTargetRegion )
                {
                  return;  // Skip this subregion
                }
              }

              string const targetFieldName = getTargetFieldName( fieldName );
              FieldSpecificationImpl::applyFieldValue< FieldSpecificationEqual >( bc, targetSet, 0.0, subRegion, targetFieldName );
            } );
          }
        }
      }
    }
  }

  return false;
}

string FieldApplicator::getTargetFieldName( string const & fieldName ) const
{
  // HydrostaticEquilibrium is defined programatically a field specification, but it isn't actually a field itself.
  // Instead, it updates the pressure field, which is the target field being modified.
  if( fieldName == "HydrostaticEquilibrium" )
  {
    return "pressure";
  }
  return fieldName;
}


void FieldApplicator::initializeSubRegionFluidState( DomainPartition & domain, ElementSubRegionBase & subRegion )
{
  GEOS_UNUSED_VAR( domain );

  // Check if this subregion has a fluid model at all
  if( !subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::fluidNamesString() ) )
  {
    return;  // No fluid model associated - nothing to initialize
  }

  // Find the CompositionalMultiphaseBase solver that targets this subregion.
  // Use m_solverName if provided, otherwise search all solvers.
  CompositionalMultiphaseBase * flowSolver = nullptr;

  Group & solversGroup = this->getGroupByPath< Group >( "/Problem/Solvers" );

  if( !m_solverName.empty() )
  {
    flowSolver = solversGroup.getGroupPointer< CompositionalMultiphaseBase >( m_solverName );
    GEOS_LOG_RANK_0_IF( flowSolver == nullptr,
                        GEOS_FMT( "FieldApplicator: solver '{}' is not a CompositionalMultiphaseBase solver; "
                                  "skipping compositional state initialization.", m_solverName ) );
  }
  else
  {
    // Search for the first CompositionalMultiphaseBase solver that has the same fluid name
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    solversGroup.forSubGroups< CompositionalMultiphaseBase >( [&]( CompositionalMultiphaseBase & solver )
    {
      if( flowSolver == nullptr )
      {
        // Check that this solver's fluid name matches the subregion's fluid
        GEOS_UNUSED_VAR( fluidName );
        flowSolver = &solver;
      }
    } );
  }

  if( flowSolver == nullptr )
  {
    // Fall back to basic fluid update only (no phaseVolFrac, relPerm, or phaseMobility)
    GEOS_LOG_RANK_0( "FieldApplicator: No CompositionalMultiphaseBase solver found. "
                     "phaseVolumeFraction, phaseMobility, and RelPerm_phaseRelPerm will not be initialized." );
    return;
  }

  // ------------------------------------------------------------------------------------
  // Step 1: Initialize porosity and permeability for the fracture (SurfaceElementSubRegion)
  //         or bulk (CellElementSubRegion).
  //
  // For SurfaceElementSubRegion (fracture elements), this computes fracPerm_permeability and
  // fracPorosity_porosity / fracPorosity_initialPorosity / fracPorosity_referencePorosity
  // from the hydraulic aperture (defaultAperture) and the reference constitutive values.
  // ------------------------------------------------------------------------------------
  if( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::solidNamesString() ) )
  {
    string const & solidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );

    if( subRegion.hasGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() ) )
    {
      dataRepository::Group & constitutiveModels =
        subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      if( constitutiveModels.hasGroup( solidName ) )
      {
        constitutive::CoupledSolidBase & porousSolid =
          constitutiveModels.getGroup< constitutive::CoupledSolidBase >( solidName );

        // Dispatch the correct updatePorosityAndPermeability overload based on subregion type
        SurfaceElementSubRegion * surfSubRegion = dynamic_cast< SurfaceElementSubRegion * >( &subRegion );
        CellElementSubRegion * cellSubRegion = dynamic_cast< CellElementSubRegion * >( &subRegion );

        if( surfSubRegion != nullptr )
        {
          // Fracture: initialize porosity/permeability from aperture (fracPerm_permeability,
          // fracPorosity_porosity, fracPorosity_initialPorosity, fracPorosity_referencePorosity)
          flowSolver->updatePorosityAndPermeability( *surfSubRegion );
        }
        else if( cellSubRegion != nullptr )
        {
          flowSolver->updatePorosityAndPermeability( *cellSubRegion );
        }

        // Initialize and save state so that _n fields (porosity_n) are correct
        porousSolid.initializeState();
        porousSolid.saveConvergedState();

        // Run update again after initializeState to ensure consistency
        if( surfSubRegion != nullptr )
        {
          flowSolver->updatePorosityAndPermeability( *surfSubRegion );
        }
        else if( cellSubRegion != nullptr )
        {
          flowSolver->updatePorosityAndPermeability( *cellSubRegion );
        }
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // Step 2: Check we have globalCompFraction and globalCompDensity to proceed with fluid init.
  // ------------------------------------------------------------------------------------
  if( !subRegion.hasWrapper( "globalCompFraction" ) || !subRegion.hasWrapper( "globalCompDensity" ) )
  {
    return;  // Not a compositional flow subregion
  }

  // ------------------------------------------------------------------------------------
  // Step 3: Initialize component densities (globalCompDensity) from component fractions.
  //         This is done via the fluid kernel wrapper (same as initializeFluidState does).
  // ------------------------------------------------------------------------------------
  {
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

    if( subRegion.hasGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() ) )
    {
      dataRepository::Group & constitutiveModels =
        subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      if( constitutiveModels.hasGroup( fluidName ) )
      {
        constitutive::MultiFluidBase & fluid =
          constitutiveModels.getGroup< constitutive::MultiFluidBase >( fluidName );

        arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
        arrayView2d< real64 const, compflow::USD_COMP > const compFrac = subRegion.getField< fields::flow::globalCompFraction >();
        arrayView2d< real64, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();

        integer const numComp = compFrac.size( 1 );

        // Update fluid model to get total density
        constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
        {
          using FluidType = TYPEOFREF( castedFluid );
          typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

          forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
          {
            fluidWrapper.update( ei, 0, pres[ei], temp[ei], compFrac[ei] );
          } );
        } );

        // Back-calculate component densities from total density and component fractions
        arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
          }
        } );
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // Step 4: Call the solver's updateFluidState to compute all derived quantities:
  //         - phaseVolumeFraction (fixes the uninitialized phaseVolumeFraction)
  //         - RelPerm_phaseRelPerm (via updateRelPermModel inside updateFluidState)
  //         - phaseMobility (via updatePhaseMobility inside updateFluidState)
  // ------------------------------------------------------------------------------------
  flowSolver->updateFluidState( subRegion );

  // ------------------------------------------------------------------------------------
  // Step 5: Save the fluid constitutive model's converged state.
  //         This saves totalDensity_n, phaseDensity_n, phaseCompFraction_n, etc.
  //         Without this, the accumulation term residual at the first Newton iteration
  //         will be non-zero (phaseDensity_n=0 leads to wrong compAmount_n).
  // ------------------------------------------------------------------------------------
  {
    string const & fluidName2 = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    if( subRegion.hasGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() ) )
    {
      dataRepository::Group & constModels2 =
        subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
      if( constModels2.hasGroup( fluidName2 ) )
      {
        constitutive::MultiFluidBase const & fluid2 =
          constModels2.getGroup< constitutive::MultiFluidBase >( fluidName2 );
        fluid2.saveConvergedState();
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // Step 6: Save the solver-level converged state so that all _n fields are correctly
  //         initialized (pressure_n, temperature_n, compDens_n, compAmount_n).
  //         This is crucial so that the first time step's accumulation term is correct.
  // ------------------------------------------------------------------------------------
  flowSolver->saveConvergedState( subRegion );

  GEOS_LOG_RANK_0( GEOS_FMT( "FieldApplicator: Fully initialized fluid state (phaseVolumeFraction, "
                             "phaseMobility, RelPerm_phaseRelPerm, fracPerm_permeability, "
                             "fracPorosity_porosity) for subregion '{}'.", subRegion.getName() ) );
}


REGISTER_CATALOG_ENTRY( TaskBase, FieldApplicator, string const &, Group * const )


} // namespace geos
