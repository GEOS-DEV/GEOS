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

  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables" );

  registerWrapper( viewKeyStruct::surfacePhaseMassDensitiesString(), &m_surfacePhaseMassDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of surface mass densities for each phase" );

  registerWrapper( viewKeyStruct::formationVolumeFactorTableNamesString(), &m_formationVolFactorTableNames ).
    setSizedFromParent( 0 );
  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setSizedFromParent( 0 );;
  registerWrapper( viewKeyStruct::waterRefPressureString(), &m_waterRefPressure );
  registerWrapper( viewKeyStruct::waterFormationVolumeFactorString(), &m_waterFormationVolFactor );
  registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility );
  registerWrapper( viewKeyStruct::waterViscosityString(), &m_waterViscosity );
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
    GEOSX_ERROR_IF( it == phaseDict.end(), "DeadOilFluid: phase not supported: " << m_phaseNames[ip] );
    integer const phaseIndex = it->second;
    GEOSX_ERROR_IF( phaseIndex >= PhaseType::MAX_NUM_PHASES, "DeadOilFluid: invalid phase index " << phaseIndex );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = LvArray::integerConversion< integer >( ip );
  }

  GEOSX_ERROR_IF( m_tableFiles.size() != numFluidPhases(),
                  "The number of table files must be equal to the number of phases" );
  for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
  {
    string const filename = m_tableFiles[ip];

    array1d< array1d< real64 > > pressureCoords;
    pressureCoords.resize( 1 );
    array1d< real64 > formationVolFactor;
    array1d< real64 > viscosity;

    std::ifstream is( filename );
    constexpr std::streamsize buf_size = 256;
    char buf[buf_size];

    localIndex counter = 0;
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

      std::cout << "strs.size() = " << strs.size() << std::endl;
      for( localIndex i = 0; i < strs.size(); ++i )
      {
        std::cout << strs[i] << std::endl;
      }

      GEOSX_ERROR_IF( (m_phaseTypes[ip] != PhaseType::WATER)
                      && strs.size() != 3,
                      "Three columns (pressure, formation volume factor, and viscosity) are expected for oil and gas" );
      GEOSX_ERROR_IF( (m_phaseTypes[ip] == PhaseType::WATER)
                      && strs.size() != 4,
                      "Four columns (pressure, formation volume factor, compressibility, and viscosity) are expected for water" );
      GEOSX_ERROR_IF( (m_phaseTypes[ip] == PhaseType::WATER)
                      && counter > 1,
                      "The water table must contain only one line" );

      if( m_phaseTypes[ip] == PhaseType::WATER )
      {
        m_waterRefPressure = stod( strs[0] );
        m_waterFormationVolFactor = stod( strs[1] );
        m_waterCompressibility = stod( strs[2] );
        m_waterViscosity = stod( strs[3] );

        GEOSX_ERROR_IF( m_waterRefPressure <= 0.0,
                        "DeadOilFluid: a strictly positive value must be provided for: "
                        << viewKeyStruct::waterRefPressureString() );
        GEOSX_ERROR_IF( m_waterFormationVolFactor <= 0.0,
                        "DeadOilFluid: a strictly positive value must be provided for: "
                        << viewKeyStruct::waterFormationVolumeFactorString() );
        GEOSX_ERROR_IF( m_waterCompressibility <= 0.0,
                        "DeadOilFluid: a strictly positive value must be provided for: "
                        << viewKeyStruct::waterCompressibilityString() );
        GEOSX_ERROR_IF( m_waterViscosity <= 0.0,
                        "DeadOilFluid: a strictly positive value must be provided for: "
                        << viewKeyStruct::waterViscosityString() );
      }
      else
      {
        pressureCoords[0].emplace_back( stod( strs[0] ) );
        formationVolFactor.emplace_back( stod( strs[1] ) );
        viscosity.emplace_back( stod( strs[2] ) );
      }

      counter++;
    }

    if( m_phaseTypes[ip] != PhaseType::WATER )
    {
      m_hydrocarbonPhaseOrder.emplace_back( LvArray::integerConversion< integer >( ip ) );
      string const formationVolFactorTableName = (m_phaseTypes[ip] == PhaseType::OIL) ? "PVDO_Bo" : "PVDG_Bg";
      string const viscosityTableName = (m_phaseTypes[ip] == PhaseType::OIL) ? "PVDO_visco" : "PVDG_viscg";
      m_formationVolFactorTableNames.emplace_back( formationVolFactorTableName );
      m_viscosityTableNames.emplace_back( viscosityTableName );

      GEOSX_ERROR_IF( pressureCoords[0].size() <= 1,
                      "DeadOilFluid: The oil and gas PVT tables must contain at least 2 values" );

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
    is.close();
  }

  // check the size of the additional parameters
  GEOSX_ERROR_IF( m_surfacePhaseMassDensity.size() != m_phaseNames.size(),
                  "DeadOilFluid: the number of surfacePhaseMassDensities is inconsistent with the phase names" );
  GEOSX_ERROR_IF( m_componentMolarWeight.size() != m_phaseNames.size(),
                  "DeadOilFluid: the number of componentMolarWeights is inconsistent with the phase names" );
}

void DeadOilFluid::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createAllKernelWrappers();
}

void DeadOilFluid::createAllKernelWrappers()
{
  FunctionManager const & functionManager = getGlobalState().getFunctionManager();

  GEOSX_ERROR_IF( m_hydrocarbonPhaseOrder.size() != 1 && m_hydrocarbonPhaseOrder.size() != 2,
                  "DeadOilFluid: the number of hydrocarbon phases should be equal to 1 (oil) or 2 (oil+gas)" );

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

  for( localIndex i = 2; i < property.size(); ++i )
  {
    GEOSX_ERROR_IF( (property[i] - property[i-1]) * (property[i-1] - property[i-2]) < 0,
                    "DeadOilFluid: the values in each PVT table must monotone" );
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
