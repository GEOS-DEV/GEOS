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
 * @file WellManager.cpp
 */

#include "WellManager.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

WellManager::WellManager( string const & name,
                          Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_useMass( false ),
  m_useTotalMassEquation( 1 ),
  m_isThermal( 0 ),
  m_isCompositional( true )
{
  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );


  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use total mass equation" );

}
Group * WellManager::createChild( string const & childKey, string const & childName )
{
  static std::set< string > const childTypes = {
    keys::compositionalMultiphaseWell,
    keys::singlePhaseWell,
    PhysicsSolverBase::groupKeyStruct::linearSolverParametersString(),
    PhysicsSolverBase::groupKeyStruct::nonlinearSolverParametersString(),
  };
  GEOS_ERROR_IF( childTypes.count( childKey ) == 0,
                 CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ),
                 getDataContext() );
  if( childKey == keys::compositionalMultiphaseWell )
  {
    setCompositional( true );
    return &registerGroup< CompositionalMultiphaseWell >( childName );
  }
  else if( childKey == keys::singlePhaseWell )
  {
    setCompositional( false );
    return &registerGroup< SinglePhaseWell >( childName );
  }
  else
  {
    PhysicsSolverBase::createChild( childKey, childName );
    return nullptr;
  }
}

void WellManager::expandObjectCatalogs()
{
  createChild( keys::compositionalMultiphaseWell, keys::compositionalMultiphaseWell );
  createChild( keys::singlePhaseWell, keys::singlePhaseWell );
}

void WellManager::registerDataOnMesh( Group & meshBodies )
{
  if( isCompositional() )
  {

    forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                      MeshLevel & mesh,
                                                      string_array const & regionNames )
    {

      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )
      {
        CompositionalMultiphaseWell & well =  getCompositionalMultiphaseWell( subRegion );
        well.registerWellDataOnMesh( subRegion );
        m_numFluidPhases = well.numFluidPhases();
        m_numFluidComponents = well.numFluidComponents();
      } );
    } );
    // 1. Set key dimensions of the problem
    // Empty check needed to avoid errors when running in schema generation mode.

    // 1 pressure + NC compositions + 1 connectionRate + temp if thermal
    m_numDofPerWellElement = isThermal() ? m_numFluidComponents + 3 : m_numFluidComponents + 2;
    // 1 pressure + NC compositions + temp if thermal
    m_numDofPerResElement = isThermal() ? m_numFluidComponents + 2 : m_numFluidComponents + 1;

  }
  else
  {
    // Single phase registration can be added here in the future
  }
#if 0
  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     WellElementSubRegion & subRegion )
    {
      if( m_referenceFluidModelName.empty() )
      {
        m_referenceFluidModelName = getConstitutiveName< MultiFluidBase >( subRegion );
      }
    } );
  } );

  // 1. Set key dimensions of the problem
  // Empty check needed to avoid errors when running in schema generation mode.
  if( !m_referenceFluidModelName.empty() )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );
    m_numPhases = fluid0.numFluidPhases();
    m_numComponents = fluid0.numFluidComponents();
  }
  // 1 pressure + NC compositions + 1 connectionRate + temp if thermal
  m_numDofPerWellElement = isThermal() ? m_numComponents + 3 : m_numComponents + 2;
  // 1 pressure + NC compositions + temp if thermal
  m_numDofPerResElement = isThermal() ? m_numComponents + 2 : m_numComponents + 1;

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerField< well::globalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< well::globalCompDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerField< well::mixtureConnectionRate >( getName() );
      subRegion.registerField< well::mixtureConnectionRate_n >( getName() );

      subRegion.registerField< well::globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< well::dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerField< well::phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< well::dPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 ); // dP, dT, dC

      subRegion.registerField< well::totalMassDensity >( getName() );
      subRegion.registerField< well::dTotalMassDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents +2 ); // dP, dT, dC

      subRegion.registerField< well::phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< well::pressureScalingFactor >( getName() );
      subRegion.registerField< well::temperatureScalingFactor >( getName() );
      subRegion.registerField< well::globalCompDensityScalingFactor >( getName() );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerField< well::compPerforationRate >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      perforationData.registerField< well::dCompPerforationRate >( getName() ).
        reference().resizeDimension< 1, 2, 3 >( 2, m_numComponents, m_numComponents+ 2 );
      if( fluid.isThermal() )
      {
        perforationData.registerField< well::energyPerforationFlux >( getName() );
        perforationData.registerField< well::dEnergyPerforationFlux >( getName() ).
          reference().resizeDimension< 1, 2 >( 2, m_numComponents+2 );
      }

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );

      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHPString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents + 2 );   // dP, dT, dC

      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );

      wellControls.registerWrapper< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0, 1 >( m_numPhases, m_numComponents + 3 );   // dP, dT, dC, dQ

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents + 3 );   // dP, dT, dC dQ

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentMassRateString() );

      // write rates output header
      // the rank that owns the reference well element is responsible
      if( m_writeCSV > 0 && subRegion.isLocallyOwned() )
      {
        string const fileName = GEOS_FMT( "{}/{}.csv", m_ratesOutputDir, wellControls.getName() );
        string const massUnit = m_useMass ? "kg" : "mol";
        integer const useSurfaceConditions = wellControls.useSurfaceConditions();
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";
        integer const numPhase = m_numPhases;
        integer const numComp = m_numComponents;
        // format: time,bhp,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
        makeDirsForPath( m_ratesOutputDir );
        GEOS_LOG( GEOS_FMT( "{}: Rates CSV generated at {}", getName(), fileName ) );
        std::ofstream outputFile( fileName );
        outputFile << "Time [s],dt[s],BHP [Pa],Total rate [" << massUnit << "/s],Total " << conditionKey << " volumetric rate [" << unitKey << "m3/s]";
        for( integer ip = 0; ip < numPhase; ++ip )
        {
          outputFile << ",Phase" << ip << " " << conditionKey << " volumetric rate [" << unitKey << "m3/s]";
        }
        for( integer ic = 0; ic < numComp; ++ic )
        {
          outputFile << ",Component" << ic << " rate [" << massUnit << "/s]";
        }
        outputFile << std::endl;
        outputFile.close();
      }
    } );
  } );
#endif
}

WellSolverBase & WellManager::getWell( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellSolverBase >( subRegion.getWellControlsName());
}

void WellManager::implicitStepSetup( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain )
{

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      WellSolverBase & well = getWell( subRegion );
      well.implicitStepSetup( time_n, dt, domain );
    } )
    ;
  } );
}

localIndex WellManager::numDofPerWellElement() const
{
  return m_numDofPerWellElement;
}

localIndex WellManager::numDofPerResElement() const
{
  return m_numDofPerResElement;
}
integer WellManager::isThermal() const
{
  return m_isThermal;
}

string WellManager::wellElementDofName() const
{
  return viewKeyStruct::dofFieldString();
}

string WellManager::resElementDofName() const
{
  return "reservoirCouplingVars";
}

localIndex WellManager::numFluidComponents() const
{
  return m_numFluidComponents;
}

localIndex WellManager::numFluidPhases() const
{
  return m_numFluidPhases;
}
WellControls & WellManager::getWellControls( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

WellControls const & WellManager::getWellControls( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

CompositionalMultiphaseWell & WellManager::getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion )
{
  return this->getGroup< CompositionalMultiphaseWell >( subRegion.getWellControlsName());
}

CompositionalMultiphaseWell const & WellManager::getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< CompositionalMultiphaseWell >( subRegion.getWellControlsName());
}
void WellManager::initializePostSubGroups()
{
#if 0
  GEOS_MARK_FUNCTION;
  // Validate constitutive models
  if( isCompositional() )
  {
    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
    constitutive::ConstitutiveManager const & cm = domain.getConstitutiveManager();
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    string const referenceFluidName = flowSolver.referenceFluidModelName();
    constitutive::MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< constitutive::MultiFluidBase >( m_referenceFluidModelName );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel const & mesh,
                                                                 string_array const & regionNames )
    {

      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion const & subRegion )
      {
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        constitutive::MultiFluidBase const & fluid = getConstitutiveModel< constitutive::MultiFluidBase >( subRegion, fluidName );
        WellControls const & wellControls = getWellControls( subRegion );
        wellControls.validateFluidModel( fluid, referenceFluid );
      } );

    } );
  }
  else
  {
    // Single phase validation can be added here in the future
  }
#endif
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, WellManager, string const &, Group * const )
} // namespace geos
