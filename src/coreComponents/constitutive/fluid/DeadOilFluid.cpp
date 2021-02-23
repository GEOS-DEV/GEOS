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

  registerWrapper( viewKeyStruct::formationVolumeFactorTableNamesString(), &m_formationVolFactorTableNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of the formation volume factor tables (one per hydrocarbon phase, in the order provided in \"phaseNames\"). \n"
                    "For instance, if \"oil\" is before \"gas\" in \"phaseNames\", the table order should be: oilTableName, gasTableName" );

  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of the viscosity tables (one per hydrocarbon phase, in the order provided in \"phaseNames\") \n"
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
    if( phaseIndex == PhaseType::OIL || phaseIndex == PhaseType::GAS )
    {
      m_hydrocarbonPhaseOrder.emplace_back( LvArray::integerConversion< integer >( ip ) );
    }
  }

  // check the water phase parameters
  integer const ipWater = m_phaseOrder[PhaseType::WATER];
  integer const ipGas = m_phaseOrder[PhaseType::GAS];
  if( ipWater >= 0 ) // if water is present
  {
    GEOSX_ERROR_IF( m_waterRefPressure <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: " << viewKeyStruct::waterRefPressureString() );
    GEOSX_ERROR_IF( m_waterFormationVolFactor <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: " << viewKeyStruct::waterFormationVolumeFactorString() );
    GEOSX_ERROR_IF( m_waterViscosity <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: " << viewKeyStruct::waterViscosityString() );
    GEOSX_ERROR_IF( m_waterCompressibility <= 0.0,
                    "DeadOilFluid: a strictly positive value must be provided for: " << viewKeyStruct::waterCompressibilityString() );
  }
  else
  {
    GEOSX_ERROR_IF( m_waterRefPressure > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed " << viewKeyStruct::waterRefPressureString() );
    GEOSX_ERROR_IF( m_waterFormationVolFactor > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed " << viewKeyStruct::waterFormationVolumeFactorString() );
    GEOSX_ERROR_IF( m_waterViscosity > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed " << viewKeyStruct::waterViscosityString() );
    GEOSX_ERROR_IF( m_waterCompressibility > 0.0,
                    "DeadOilFluid: if water is absent, this keyword is not needed " << viewKeyStruct::waterCompressibilityString() );
  }

  // check the table names
  localIndex const numExpectedTables = (ipGas >= 0) ? 2 : 1;
  GEOSX_ERROR_IF( m_formationVolFactorTableNames.size() != numExpectedTables,
                  "DeadOilFluid: one formation volume factor table must be provided for each hydrocarbon phase" );
  GEOSX_ERROR_IF( m_viscosityTableNames.size() != numExpectedTables,
                  "DeadOilFluid: one viscosity table must be provided for each hydrocarbon phase" );

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
      localIndex const ip = m_hydrocarbonPhaseOrder[iph];

      // grab the tables by name from the function manager
      TableFunction const & FVFTable = functionManager.getGroup< TableFunction const >( m_formationVolFactorTableNames[ip] );
      TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[ip] );
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

  GEOSX_ERROR_IF( table.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                  "DeadOilFluid: the interpolation method for the PVT tables must be linear" );
  GEOSX_ERROR_IF( property.size() <= 1,
                  "DeadOilFluid: each PVT table must contain at least 2 values" );

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
