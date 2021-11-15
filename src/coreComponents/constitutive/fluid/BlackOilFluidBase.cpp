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

#include "BlackOilFluidBase.hpp"

#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

BlackOilFluidBase::BlackOilFluidBase( string const & name,
                                      Group * const parent )
  :
  MultiFluidBase( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::surfacePhaseMassDensitiesString(), &m_surfacePhaseMassDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of surface mass densities for each phase" );

  // 1) First option: specify PVT tables from one file per phase, read the files line by line, and populate the internal TableFunctions

  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables (one per phase)" );

  // 2) Second option: specify TableFunction names for oil and gas,
  //                   and then specify water reference pressure, water Bw, water compressibility and water viscosity

  registerWrapper( viewKeyStruct::formationVolumeFactorTableNamesString(), &m_formationVolFactorTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of formation volume factor TableFunction names from the Functions block. \n"
                    "The user must provide one TableFunction per hydrocarbon phase, in the order provided in \"phaseNames\". \n"
                    "For instance, if \"oil\" is before \"gas\" in \"phaseNames\", the table order should be: oilTableName, gasTableName" );
  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of viscosity TableFunction names from the Functions block. \n"
                    "The user must provide one TableFunction per hydrocarbon phase, in the order provided in \"phaseNames\". \n"
                    "For instance, if \"oil\" is before \"gas\" in \"phaseNames\", the table order should be: oilTableName, gasTableName" );
  registerWrapper( viewKeyStruct::waterRefPressureString(), &m_waterParams.referencePressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water reference pressure" );
  registerWrapper( viewKeyStruct::waterFormationVolumeFactorString(), &m_waterParams.formationVolFactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water formation volume factor" );
  registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterParams.compressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water compressibility" );
  registerWrapper( viewKeyStruct::waterViscosityString(), &m_waterParams.viscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water viscosity" );

  // Register extra wrappers to enable auto-cloning
  registerWrapper( "phaseTypes", &m_phaseTypes )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "phaseOrder", &m_phaseOrder )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "hydrocarbonPhaseOrder", &m_hydrocarbonPhaseOrder )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "formationVolFactorTableWrappers", &m_formationVolFactorTables )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "viscosityTableWrappers", &m_viscosityTables )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
}

void BlackOilFluidBase::fillWaterData( array1d< array1d< real64 > > const & tableValues )
{
  GEOSX_THROW_IF_NE_MSG( tableValues.size(), 1,
                         getFullName() << ": the water table must contain one line and only one",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( tableValues[0].size(), 4,
                         getFullName() << ": four columns (pressure, formation volume factor, compressibility, and viscosity) are expected for water",
                         InputError );

  GEOSX_THROW_IF( m_waterParams.referencePressure > 0.0 || m_waterParams.formationVolFactor > 0.0 ||
                  m_waterParams.compressibility > 0.0 || m_waterParams.viscosity > 0.0,
                  getFullName() << ": input is redundant (user provided both water data and a water pvt file)",
                  InputError );

  m_waterParams.referencePressure = tableValues[0][0];
  m_waterParams.formationVolFactor = tableValues[0][1];
  m_waterParams.compressibility = tableValues[0][2];
  m_waterParams.viscosity = tableValues[0][3];

  validateWaterParams();
}

void BlackOilFluidBase::fillHydrocarbonData( integer const ip,
                                             array1d< array1d< real64 > > const & tableValues )
{
  array1d< array1d< real64 > > pressureCoords( 1 );
  pressureCoords[0].resize( tableValues.size() );
  array1d< real64 > formationVolFactor( tableValues.size() );
  array1d< real64 > viscosity( tableValues.size() );

  for( localIndex i = 0; i < tableValues.size(); ++i )
  {
    GEOSX_THROW_IF_NE_MSG( tableValues[i].size(), 3,
                           getFullName() << ": three columns (pressure, formation volume factor, and viscosity) are expected for oil and gas",
                           InputError );

    pressureCoords[0][i] = tableValues[i][0];
    formationVolFactor[i] = tableValues[i][1];
    viscosity[i] = tableValues[i][2];
  }

  m_hydrocarbonPhaseOrder.emplace_back( ip );
  string const formationVolFactorTableName = getName() + (m_phaseTypes[ip] == PhaseType::OIL ? "_PVDO_Bo" : "_PVDG_Bg");
  string const viscosityTableName = getName() +  (m_phaseTypes[ip] == PhaseType::OIL ? "_PVDO_visco" : "_PVDG_viscg");
  m_formationVolFactorTableNames.emplace_back( formationVolFactorTableName );
  m_viscosityTableNames.emplace_back( viscosityTableName );

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & tablePVDX_B =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", formationVolFactorTableName ) );
  tablePVDX_B.setTableCoordinates( pressureCoords );
  tablePVDX_B.setTableValues( formationVolFactor );
  tablePVDX_B.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDX_visc =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", viscosityTableName ) );
  tablePVDX_visc.setTableCoordinates( pressureCoords );
  tablePVDX_visc.setTableValues( viscosity );
  tablePVDX_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

integer BlackOilFluidBase::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

void BlackOilFluidBase::postProcessInput()
{
  m_componentNames = m_phaseNames;

  MultiFluidBase::postProcessInput();

  m_phaseTypes.resize( numFluidPhases() );
  m_phaseOrder.resizeDefault( MAX_NUM_PHASES, -1 );

  auto const toPhaseType = [&]( string const & lookup )
  {
    static unordered_map< string, integer > const phaseDict =
    {
      { "gas", BlackOilFluidBase::PhaseType::GAS },
      { "oil", BlackOilFluidBase::PhaseType::OIL },
      { "water", BlackOilFluidBase::PhaseType::WATER }
    };
    return findOption( phaseDict, lookup, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    m_phaseTypes[ip] = toPhaseType( m_phaseNames[ip] );
    m_phaseOrder[m_phaseTypes[ip]] = ip;
  }

  auto const checkInputSize = [&]( auto const & array, auto const & attribute )
  {
    GEOSX_THROW_IF_NE_MSG( array.size(), m_phaseNames.size(),
                           GEOSX_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                           InputError );
  };
  checkInputSize( m_surfacePhaseMassDensity, viewKeyStruct::surfacePhaseMassDensitiesString() );
  checkInputSize( m_componentMolarWeight, viewKeyStruct::componentMolarWeightString() );

  // we make the distinction between the two input options
  if( m_tableFiles.empty() )
  {
    readInputDataFromTableFunctions();
  }
  else
  {
    readInputDataFromPVTFiles();
  }

}

void BlackOilFluidBase::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createAllKernelWrappers();
}

void BlackOilFluidBase::createAllKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  GEOSX_THROW_IF( m_hydrocarbonPhaseOrder.size() != 1 && m_hydrocarbonPhaseOrder.size() != 2,
                  GEOSX_FMT( "{}: the number of hydrocarbon phases must be 1 (oil) or 2 (oil+gas)", getFullName() ),
                  InputError );

  if( m_formationVolFactorTables.empty() && m_viscosityTables.empty() )
  {

    // loop over the hydrocarbon phases
    for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
    {
      // grab the tables by name from the function manager
      TableFunction const & fvfTable = functionManager.getGroup< TableFunction const >( m_formationVolFactorTableNames[iph] );
      TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[iph] );
      validateTable( fvfTable );
      validateTable( viscosityTable );

      // create the table wrapper for the oil and (if present) the gas phases
      m_formationVolFactorTables.emplace_back( fvfTable.createKernelWrapper() );
      m_viscosityTables.emplace_back( viscosityTable.createKernelWrapper() );
    }
  }
}

void BlackOilFluidBase::validateTable( TableFunction const & table ) const
{
  arrayView1d< real64 const > const property = table.getValues();

  GEOSX_THROW_IF_NE_MSG( table.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                         GEOSX_FMT( "{}: in table '{}' interpolation method must be linear", getFullName(), table.getName() ),
                         InputError );
  GEOSX_THROW_IF_LT_MSG( property.size(), 2,
                         GEOSX_FMT( "{}: table `{}` must contain at least two values", getFullName(), table.getName() ),
                         InputError );

  for( localIndex i = 2; i < property.size(); ++i )
  {
    GEOSX_THROW_IF( (property[i] - property[i-1]) * (property[i-1] - property[i-2]) < 0,
                    GEOSX_FMT( "{}: in table '{}' values must be monotone", getFullName(), table.getName() ),
                    InputError );
  }
}

void BlackOilFluidBase::validateWaterParams() const
{
  auto const checkPositiveValue = [&]( real64 const value, auto const & attribute )
  {
    GEOSX_THROW_IF_LE_MSG( value, 0.0,
                           GEOSX_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                           InputError );
  };
  checkPositiveValue( m_waterParams.referencePressure, viewKeyStruct::waterRefPressureString() );
  checkPositiveValue( m_waterParams.formationVolFactor, viewKeyStruct::waterFormationVolumeFactorString() );
  checkPositiveValue( m_waterParams.compressibility, viewKeyStruct::waterCompressibilityString() );
  checkPositiveValue( m_waterParams.viscosity, viewKeyStruct::waterViscosityString() );
}

BlackOilFluidBase::KernelWrapper::
  KernelWrapper( arrayView1d< integer const > phaseTypes,
                 arrayView1d< integer const > phaseOrder,
                 arrayView1d< integer const > hydrocarbonPhaseOrder,
                 arrayView1d< real64 const > surfacePhaseMassDensity,
                 arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                 arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                 BlackOilFluidBase::WaterParams const waterParams,
                 arrayView1d< real64 const > componentMolarWeight,
                 bool const useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_phaseTypes( std::move( phaseTypes ) ),
  m_phaseOrder( std::move( phaseOrder ) ),
  m_hydrocarbonPhaseOrder( std::move( hydrocarbonPhaseOrder ) ),
  m_surfacePhaseMassDensity( std::move( surfacePhaseMassDensity ) ),
  m_formationVolFactorTables( std::move( formationVolFactorTables ) ),
  m_viscosityTables( std::move( viscosityTables ) ),
  m_waterParams( waterParams )
{}

} // namespace constitutive

} // namespace geosx
