/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "DeadOilFluid.hpp"

#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "managers/GeosxState.hpp"
#include "managers/Functions/FunctionManager.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace stringutilities;

namespace constitutive
{

constexpr integer DeadOilFluid::PhaseType::GAS;
constexpr integer DeadOilFluid::PhaseType::OIL;
constexpr integer DeadOilFluid::PhaseType::WATER;

namespace
{

std::unordered_map< string, integer > const phaseDict =
{
  { "gas", DeadOilFluid::PhaseType::GAS   },
  { "oil", DeadOilFluid::PhaseType::OIL   },
  { "water", DeadOilFluid::PhaseType::WATER }
};

}

DeadOilFluid::DeadOilFluid( string const & name,
                            Group * const parent )
  :
  MultiFluidBase( name, parent ),
  m_waterRefPressure( 0.0 ),
  m_waterFormationVolFactor( 0.0 ),
  m_waterCompressibility( 0.0 ),
  m_waterViscosity( 0.0 )
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
  registerWrapper( viewKeyStruct::waterRefPressureString(), &m_waterRefPressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water reference pressure" );
  registerWrapper( viewKeyStruct::waterFormationVolumeFactorString(), &m_waterFormationVolFactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water formation volume factor" );
  registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water compressibility" );
  registerWrapper( viewKeyStruct::waterViscosityString(), &m_waterViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water viscosity" );
}

void DeadOilFluid::readTable( string const & filename,
                              array1d< array1d< real64 > > & values ) const
{
  std::ifstream is( filename );
  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  while( is.getline( buf, buf_size ) )
  {
    string const str( buf );
    string_array const strs = Tokenize( str, " " );

    if( strs.empty() )
    {
      continue;
    }
    if( strs[0].front() == '#' )
    {
      continue;
    }
    if( streq( strs[0], "--" ) )
    {
      continue;
    }

    array1d< real64 > rowValues;
    rowValues.resize( strs.size() );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      rowValues[i] = stod( strs[i] );
    }
    values.emplace_back( rowValues );
  }
  is.close();
}

void DeadOilFluid::fillWaterData( array1d< array1d< real64 > > const & tableValues )
{
  GEOSX_THROW_IF( tableValues.size() != 1,
                  "DeadOilFluid: the water table must contain one line and only one",
                  InputError );
  GEOSX_THROW_IF( tableValues[0].size() != 4,
                  "DeadOilFluid: four columns (pressure, formation volume factor, compressibility, and viscosity) are expected for water",
                  InputError );

  GEOSX_THROW_IF( m_waterRefPressure > 0.0 || m_waterFormationVolFactor > 0.0
                  || m_waterCompressibility > 0.0 || m_waterViscosity > 0.0,
                  "DeadOilFluid: input is redundant (user provided both water data and a water pvt file)",
                  InputError );

  m_waterRefPressure = tableValues[0][0];
  m_waterFormationVolFactor = tableValues[0][1];
  m_waterCompressibility = tableValues[0][2];
  m_waterViscosity = tableValues[0][3];

  GEOSX_THROW_IF( m_waterRefPressure <= 0.0,
                  "DeadOilFluid: a strictly positive value must be provided for: "
                  << viewKeyStruct::waterRefPressureString(),
                  InputError );
  GEOSX_THROW_IF( m_waterFormationVolFactor <= 0.0,
                  "DeadOilFluid: a strictly positive value must be provided for: "
                  << viewKeyStruct::waterFormationVolumeFactorString(),
                  InputError );
  GEOSX_THROW_IF( m_waterCompressibility <= 0.0,
                  "DeadOilFluid: a strictly positive value must be provided for: "
                  << viewKeyStruct::waterCompressibilityString(),
                  InputError );
  GEOSX_THROW_IF( m_waterViscosity <= 0.0,
                  "DeadOilFluid: a strictly positive value must be provided for: "
                  << viewKeyStruct::waterViscosityString(),
                  InputError );
}

void DeadOilFluid::fillHydrocarbonData( localIndex const ip,
                                        array1d< array1d< real64 > > const & tableValues )
{
  array1d< array1d< real64 > > pressureCoords;
  pressureCoords.resize( 1 );
  pressureCoords[0].resize( tableValues.size() );
  array1d< real64 > formationVolFactor;
  formationVolFactor.resize( tableValues.size() );
  array1d< real64 > viscosity;
  viscosity.resize( tableValues.size() );

  for( localIndex i = 0; i < tableValues.size(); ++i )
  {
    GEOSX_THROW_IF( tableValues[i].size() != 3,
                    "DeadOilFluid: three columns (pressure, formation volume factor, and viscosity) are expected for oil and gas",
                    InputError );

    pressureCoords[0][i] = tableValues[i][0];
    formationVolFactor[i] = tableValues[i][1];
    viscosity[i] = tableValues[i][2];
  }

  m_hydrocarbonPhaseOrder.emplace_back( LvArray::integerConversion< integer >( ip ) );
  string const formationVolFactorTableName = (m_phaseTypes[ip] == PhaseType::OIL) ? "PVDO_Bo" : "PVDG_Bg";
  string const viscosityTableName = (m_phaseTypes[ip] == PhaseType::OIL) ? "PVDO_visco" : "PVDG_viscg";
  m_formationVolFactorTableNames.emplace_back( formationVolFactorTableName );
  m_viscosityTableNames.emplace_back( viscosityTableName );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  TableFunction & tablePVDX_B =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", formationVolFactorTableName ) );
  tablePVDX_B.setTableCoordinates( pressureCoords );
  tablePVDX_B.setTableValues( formationVolFactor );
  tablePVDX_B.reInitializeFunction();
  tablePVDX_B.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  TableFunction & tablePVDX_visc =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", viscosityTableName ) );
  tablePVDX_visc.setTableCoordinates( pressureCoords );
  tablePVDX_visc.setTableValues( viscosity );
  tablePVDX_visc.reInitializeFunction();
  tablePVDX_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void DeadOilFluid::readInputDataFromPVTFiles()
{
  GEOSX_THROW_IF( m_tableFiles.size() != numFluidPhases(),
                  "DeadOilFluid: the number of table files must be equal to the number of phases",
                  InputError );
  GEOSX_THROW_IF( m_formationVolFactorTableNames.size() > 0.0 || m_viscosityTableNames.size() > 0.0,
                  "DeadOilFluid: input is redundant (user provided both TableFunction names and pvt files)",
                  InputError );

  array1d< array1d< real64 > > tableValues;
  for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
  {
    tableValues.clear();
    readTable( m_tableFiles[ip], tableValues );

    if( m_phaseTypes[ip] == PhaseType::WATER )
    {
      fillWaterData( tableValues );
    }
    else
    {
      fillHydrocarbonData( ip, tableValues );
    }
  }
}

void DeadOilFluid::useProvidedTableFunctions()
{
  GEOSX_THROW_IF( m_tableFiles.size() > 0,
                  "DeadOilFluid: input is redundant (user provided both TableFunction names and pvt files)",
                  InputError );

  integer const ipWater = m_phaseOrder[PhaseType::WATER];
  integer const ipGas = m_phaseOrder[PhaseType::GAS];
  if( ipWater >= 0 ) // if water is present
  {
    GEOSX_THROW_IF( m_waterRefPressure <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: "
                    << viewKeyStruct::waterRefPressureString(),
                    InputError );
    GEOSX_THROW_IF( m_waterFormationVolFactor <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: "
                    << viewKeyStruct::waterFormationVolumeFactorString(),
                    InputError );
    GEOSX_THROW_IF( m_waterCompressibility <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: "
                    << viewKeyStruct::waterCompressibilityString(),
                    InputError );
    GEOSX_THROW_IF( m_waterViscosity <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: "
                    << viewKeyStruct::waterViscosityString(),
                    InputError );
  }
  else
  {
    GEOSX_THROW_IF( m_waterRefPressure > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed "
                    << viewKeyStruct::waterRefPressureString(),
                    InputError );
    GEOSX_THROW_IF( m_waterFormationVolFactor > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed "
                    << viewKeyStruct::waterFormationVolumeFactorString(),
                    InputError );
    GEOSX_THROW_IF( m_waterViscosity > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed "
                    << viewKeyStruct::waterViscosityString(),
                    InputError );
    GEOSX_THROW_IF( m_waterCompressibility > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed "
                    << viewKeyStruct::waterCompressibilityString(),
                    InputError );
  }

  localIndex const numExpectedTables = (ipGas >= 0) ? 2 : 1;
  GEOSX_THROW_IF( m_formationVolFactorTableNames.size() != numExpectedTables,
                  "DeadOilFluid: one formation volume factor table must be provided for each hydrocarbon phase",
                  InputError );
  GEOSX_THROW_IF( m_viscosityTableNames.size() != numExpectedTables,
                  "DeadOilFluid: one viscosity table must be provided for each hydrocarbon phase",
                  InputError );

  for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
  {
    if( m_phaseTypes[ip] == PhaseType::OIL || m_phaseTypes[ip] == PhaseType::GAS )
    {
      m_hydrocarbonPhaseOrder.emplace_back( LvArray::integerConversion< integer >( ip ) );
    }
  }

  FunctionManager const & functionManager = getGlobalState().getFunctionManager();
  for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    GEOSX_THROW_IF( !functionManager.hasGroup( m_formationVolFactorTableNames[iph] ),
                    "DeadOilFluid: the formation volume factor table " << m_formationVolFactorTableNames[iph]
                                                                       << " could not be found",
                    InputError );
    GEOSX_THROW_IF( !functionManager.hasGroup( m_viscosityTableNames[iph] ),
                    "DeadOilFluid: the viscosity table " << m_viscosityTableNames[iph]
                                                         << " could not be found",
                    InputError );
  }
}

void DeadOilFluid::postProcessInput()
{
  m_componentNames = m_phaseNames;

  MultiFluidBase::postProcessInput();

  m_phaseTypes.resize( numFluidPhases() );
  m_phaseOrder.resize( PhaseType::MAX_NUM_PHASES );
  m_phaseOrder.setValues< serialPolicy >( -1 );

  for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
  {
    auto it = phaseDict.find( m_phaseNames[ip] );
    GEOSX_THROW_IF( it == phaseDict.end(),
                    "DeadOilFluid: phase not supported: " << m_phaseNames[ip],
                    InputError );
    integer const phaseIndex = it->second;
    GEOSX_THROW_IF( phaseIndex >= PhaseType::MAX_NUM_PHASES, "DeadOilFluid: invalid phase index " << phaseIndex,
                    InputError );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = LvArray::integerConversion< integer >( ip );
  }

  GEOSX_THROW_IF( m_surfacePhaseMassDensity.size() != m_phaseNames.size(),
                  "DeadOilFluid: the number of surfacePhaseMassDensities is inconsistent with the phase names",
                  InputError );
  GEOSX_THROW_IF( m_componentMolarWeight.size() != m_phaseNames.size(),
                  "DeadOilFluid: the number of componentMolarWeights is inconsistent with the phase names",
                  InputError );

  // at this point, we make the distinction between the two input options
  if( !m_tableFiles.empty() )
  {
    readInputDataFromPVTFiles();
  }
  else
  {
    useProvidedTableFunctions();
  }
}

void DeadOilFluid::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createAllKernelWrappers();
}

void DeadOilFluid::createAllKernelWrappers()
{
  FunctionManager const & functionManager = getGlobalState().getFunctionManager();

  GEOSX_THROW_IF( m_hydrocarbonPhaseOrder.size() != 1 && m_hydrocarbonPhaseOrder.size() != 2,
                  "DeadOilFluid: the number of hydrocarbon phases should be equal to 1 (oil) or 2 (oil+gas)",
                  InputError );

  if( m_formationVolFactorTables.size() == 0 && m_viscosityTables.size() == 0 )
  {

    // loop over the hydrocarbon phases
    for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
    {
      // grab the tables by name from the function manager
      TableFunction const & FVFTable = functionManager.getGroup< TableFunction const >( m_formationVolFactorTableNames[iph] );
      TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[iph] );
      validateTable( FVFTable );
      validateTable( viscosityTable );

      // create the table wrapper for the oil and (if present) the gas phases
      m_formationVolFactorTables.emplace_back( FVFTable.createKernelWrapper() );
      m_viscosityTables.emplace_back( viscosityTable.createKernelWrapper() );
    }
  }

}

void DeadOilFluid::validateTable( TableFunction const & table ) const
{
  array1d< real64 > const & property = table.getValues();

  GEOSX_THROW_IF( table.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                  "DeadOilFluid: the interpolation method in the table must be linear",
                  InputError );
  GEOSX_THROW_IF( property.size() <= 1,
                  "DeadOilFluid: each PVT table must contain at least 2 values",
                  InputError );

  for( localIndex i = 2; i < property.size(); ++i )
  {
    GEOSX_THROW_IF( (property[i] - property[i-1]) * (property[i-1] - property[i-2]) < 0,
                    "DeadOilFluid: the values in each PVT table must monotone",
                    InputError );
  }
}

std::unique_ptr< ConstitutiveBase >
DeadOilFluid::deliverClone( string const & name,
                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  DeadOilFluid & model = dynamicCast< DeadOilFluid & >( *clone );

  model.m_phaseTypes = m_phaseTypes;
  model.m_phaseOrder = m_phaseOrder;
  model.m_hydrocarbonPhaseOrder = m_hydrocarbonPhaseOrder;

  model.createAllKernelWrappers();

  return clone;
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, DeadOilFluid, string const &, Group * const )
} //namespace constitutive

} //namespace geosx
