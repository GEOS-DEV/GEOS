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
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace stringutilities;
namespace constitutive
{

DeadOilFluid::DeadOilFluid( string const & name,
                            Group * const parent )
  :
  BlackOilFluidBase( name, parent )
{}

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
    PVTProps::BlackOilTables::readTable( m_tableFiles[ip], 3, tableValues );

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

  FunctionManager const & functionManager = FunctionManager::getInstance();
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
