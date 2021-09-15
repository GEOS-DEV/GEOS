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
  GEOSX_THROW_IF_NE_MSG( m_tableFiles.size(), numFluidPhases(),
                         getFullName() << ": the number of table files must be equal to the number of phases",
                         InputError );
  GEOSX_THROW_IF( m_formationVolFactorTableNames.size() > 0.0 || m_viscosityTableNames.size() > 0.0,
                  getFullName() << ": input is redundant (user provided both TableFunction names and pvt files)",
                  InputError );

  array1d< array1d< real64 > > tableValues;
  for( integer ip = 0; ip < numFluidPhases(); ++ip )
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

void DeadOilFluid::readInputDataFromTableFunctions()
{
  GEOSX_THROW_IF( !m_tableFiles.empty(),
                  getFullName() << ": input is redundant (user provided both TableFunction names and pvt files)",
                  InputError );

  integer const ipWater = m_phaseOrder[PhaseType::WATER];
  integer const ipGas = m_phaseOrder[PhaseType::GAS];
  if( ipWater >= 0 ) // if water is present
  {
    validateWaterParams();
  }
  else
  {
    GEOSX_THROW_IF( m_waterParams.referencePressure > 0.0,
                    getFullName() << ": if water is absent, this keyword is not needed " << viewKeyStruct::waterRefPressureString(),
                    InputError );
    GEOSX_THROW_IF( m_waterParams.formationVolFactor > 0.0,
                    getFullName() << ": if water is absent, this keyword is not needed " << viewKeyStruct::waterFormationVolumeFactorString(),
                    InputError );
    GEOSX_THROW_IF( m_waterParams.viscosity > 0.0,
                    getFullName() << ": if water is absent, this keyword is not needed " << viewKeyStruct::waterViscosityString(),
                    InputError );
    GEOSX_THROW_IF( m_waterParams.compressibility > 0.0,
                    getFullName() << ": if water is absent, this keyword is not needed " << viewKeyStruct::waterCompressibilityString(),
                    InputError );
  }

  integer const numExpectedTables = (ipGas >= 0) ? 2 : 1;
  GEOSX_THROW_IF_NE_MSG( m_formationVolFactorTableNames.size(), numExpectedTables,
                         getFullName() << ": one formation volume factor table must be provided for each hydrocarbon phase",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( m_viscosityTableNames.size(), numExpectedTables,
                         getFullName() << ": one viscosity table must be provided for each hydrocarbon phase",
                         InputError );

  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    if( m_phaseTypes[ip] == PhaseType::OIL || m_phaseTypes[ip] == PhaseType::GAS )
    {
      m_hydrocarbonPhaseOrder.emplace_back( ip );
    }
  }

  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    GEOSX_THROW_IF( !functionManager.hasGroup( m_formationVolFactorTableNames[iph] ),
                    getFullName() << ": the formation volume factor table " << m_formationVolFactorTableNames[iph] << " could not be found",
                    InputError );
    GEOSX_THROW_IF( !functionManager.hasGroup( m_viscosityTableNames[iph] ),
                    getFullName() << ": the viscosity table " << m_viscosityTableNames[iph] << " could not be found",
                    InputError );
  }
}

DeadOilFluid::KernelWrapper::
  KernelWrapper( arrayView1d< integer const > phaseTypes,
                 arrayView1d< integer const > phaseOrder,
                 arrayView1d< integer const > hydrocarbonPhaseOrder,
                 arrayView1d< real64 const > surfacePhaseMassDensity,
                 arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                 arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                 BlackOilFluidBase::WaterParams const waterParams,
                 arrayView1d< real64 const > componentMolarWeight,
                 bool useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : BlackOilFluidBase::KernelWrapper( std::move( phaseTypes ),
                                      std::move( phaseOrder ),
                                      std::move( hydrocarbonPhaseOrder ),
                                      std::move( surfacePhaseMassDensity ),
                                      std::move( formationVolFactorTables ),
                                      std::move( viscosityTables ),
                                      waterParams,
                                      std::move( componentMolarWeight ),
                                      useMass,
                                      std::move( phaseFraction ),
                                      std::move( phaseDensity ),
                                      std::move( phaseMassDensity ),
                                      std::move( phaseViscosity ),
                                      std::move( phaseCompFraction ),
                                      std::move( totalDensity ) )
{}

DeadOilFluid::KernelWrapper
DeadOilFluid::createKernelWrapper()
{
  return KernelWrapper( m_phaseTypes,
                        m_phaseOrder,
                        m_hydrocarbonPhaseOrder,
                        m_surfacePhaseMassDensity,
                        m_formationVolFactorTables,
                        m_viscosityTables,
                        m_waterParams,
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DeadOilFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
