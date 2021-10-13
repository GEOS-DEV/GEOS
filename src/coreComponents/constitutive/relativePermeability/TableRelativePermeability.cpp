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

/**
 * @file TableRelativePermeability.cpp
 */

#include "TableRelativePermeability.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableRelativePermeability::TableRelativePermeability( std::string const & name,
                                                      Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::wettingNonWettingRelPermTableNamesString(), &m_wettingNonWettingRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (wetting phase, non-wetting phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, nonWettingPhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateRelPermTableNamesString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::wettingIntermediateRelPermTableNamesString(), &m_wettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (wetting phase, intermediate phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::nonWettingIntermediateRelPermTableNamesString(), &m_nonWettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (non-wetting phase, intermediate phase)\n"
                    "The expected format is \"{ nonWettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( "relPermWrappers", &m_relPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TableRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  localIndex const numPhases = m_phaseNames.size();
  GEOSX_THROW_IF( numPhases != 2 && numPhases != 3,
                  GEOSX_FMT( "{}: the expected number of fluid phases is either two, or three",
                             getFullName() ),
                  InputError );

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( m_wettingNonWettingRelPermTableNames.empty(),
                    GEOSX_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the relative permeability tables for the pair (wetting phase, non-wetting phase)",
                               getFullName(),
                               viewKeyStruct::wettingNonWettingRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_wettingNonWettingRelPermTableNames.size() != 2,
                    GEOSX_FMT(
                      "{}: for a two-phase flow simulation, we must use {} to specify exactly two names: first the name of the wetting phase relperm table, second the name on the non-wetting phase relperm table",
                      getFullName(),
                      viewKeyStruct::wettingNonWettingRelPermTableNamesString() ),
                    InputError );

  }
  else if( numPhases == 3 )
  {
    GEOSX_THROW_IF( m_wettingIntermediateRelPermTableNames.empty() || m_nonWettingIntermediateRelPermTableNames.empty(),
                    GEOSX_FMT(
                      "{}: for a three-phase flow simulation, we must use {} to specify the relative permeability tables for the pair (wetting phase, intermediate phase), and {} to specify the relative permeability tables for the pair (non-wetting phase, intermediate phase)",
                      getFullName(),
                      viewKeyStruct::wettingIntermediateRelPermTableNamesString(),
                      viewKeyStruct::nonWettingIntermediateRelPermTableNamesString()  ),
                    InputError );

    GEOSX_THROW_IF( m_wettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT(
                      "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: first the name of the wetting phase relperm table, second the name on the intermediate phase relperm table",
                      getFullName(),
                      viewKeyStruct::wettingIntermediateRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_nonWettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT(
                      "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: first the name of the non-wetting phase relperm table, second the name on the intermediate phase relperm table",
                      getFullName(),
                      viewKeyStruct::nonWettingIntermediateRelPermTableNamesString() ),
                    InputError );
  }
}

void TableRelativePermeability::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();

  localIndex const numPhases = m_phaseNames.size();
  m_phaseMinVolumeFraction.resize( MAX_NUM_PHASES );

  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( numPhases == 2 )
  {
    for( integer ip = 0; ip < m_wettingNonWettingRelPermTableNames.size(); ++ip )
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_wettingNonWettingRelPermTableNames[ip] ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_wettingNonWettingRelPermTableNames[ip] ),
                      InputError );
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingRelPermTableNames[ip] );
      real64 const minVolPhaseFrac = validateRelativePermeabilityTable( relPermTable );
      if( ip == 0 ) // wetting phase is either water, or oil (for two-phase oil-gas systems)
      {
        localIndex const ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
        m_phaseMinVolumeFraction[ipWetting] = minVolPhaseFrac;
      }
      else if( ip == 1 ) // non-wetting phase is either oil (for two-phase oil-water systems), or gas
      {
        localIndex const ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
        m_phaseMinVolumeFraction[ipNonWetting] = minVolPhaseFrac;
      }
    }
  }
  else if( numPhases == 3 )
  {
    for( integer ip = 0; ip < m_wettingIntermediateRelPermTableNames.size(); ++ip )
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_wettingIntermediateRelPermTableNames[ip] ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_wettingIntermediateRelPermTableNames[ip] ),
                      InputError );
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_wettingIntermediateRelPermTableNames[ip] );
      real64 const minVolPhaseFrac = validateRelativePermeabilityTable( relPermTable );

      if( ip == 0 ) // wetting phase is water
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::WATER]] = minVolPhaseFrac;
      }
      else if( ip == 1 ) // intermediate phase is oil
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::OIL]] = minVolPhaseFrac;
      }
    }
    for( integer ip = 0; ip < m_nonWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_nonWettingIntermediateRelPermTableNames[ip] ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_nonWettingIntermediateRelPermTableNames[ip] ),
                      InputError );
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateRelPermTableNames[ip] );
      real64 const minVolPhaseFrac = validateRelativePermeabilityTable( relPermTable );

      if( ip == 0 ) // non-wetting phase is gas
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::GAS]] = minVolPhaseFrac;
      }
      else if( ip == 1 ) // intermediate phase is oil
      {
        m_phaseMinVolumeFraction[m_phaseOrder[PhaseType::OIL]] = minVolPhaseFrac;
      }
    }
  }
}

void TableRelativePermeability::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  localIndex const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_relPermKernelWrappers.clear();
  if( numPhases == 2 )
  {
    for( integer ip = 0; ip < m_wettingNonWettingRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingRelPermTableNames[ip] );
      m_relPermKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }
  }
  else if( numPhases == 3 )
  {
    for( integer ip = 0; ip < m_wettingIntermediateRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_wettingIntermediateRelPermTableNames[ip] );
      m_relPermKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }
    for( integer ip = 0; ip < m_nonWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateRelPermTableNames[ip] );
      m_relPermKernelWrappers.emplace_back( relPermTable.createKernelWrapper() );
    }
  }
}

real64 TableRelativePermeability::validateRelativePermeabilityTable( TableFunction const & relPermTable ) const
{
  ArrayOfArraysView< real64 const > coords = relPermTable.getCoordinates();

  GEOSX_THROW_IF_NE_MSG( relPermTable.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                         GEOSX_FMT( "{}: in table '{}' interpolation method must be linear", getFullName(), relPermTable.getName() ),
                         InputError );
  GEOSX_THROW_IF_NE_MSG( relPermTable.numDimensions(), 1,
                         GEOSX_FMT( "{}: table '{}' must have a single independent coordinate", getFullName(), relPermTable.getName() ),
                         InputError );
  GEOSX_THROW_IF_LT_MSG( coords.sizeOfArray( 0 ), 2,
                         GEOSX_FMT( "{}: table `{}` must contain at least two values", getFullName(), relPermTable.getName() ),
                         InputError );

  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const relPerm = relPermTable.getValues();
  real64 minVolFraction = phaseVolFrac[0];

  // note that the TableFunction class has already checked that coords.sizeOfArray( 0 ) == relPerm.size()
  GEOSX_THROW_IF( !isZero( relPerm[0] ),
                  GEOSX_FMT( "{}: in table '{}' the first value must be equal to 0", getFullName(), relPermTable.getName() ),
                  InputError );
  for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
  {
    // check phase volume fraction
    GEOSX_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                    GEOSX_FMT( "{}: in table '{}' values must be between 0 and 1", getFullName(), relPermTable.getName() ),
                    InputError );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check phase relative permeability
    GEOSX_THROW_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-10,
                    GEOSX_FMT( "{}: in table '{}' values must be strictly increasing", getFullName(), relPermTable.getName() ),
                    InputError );

    if( isZero( relPerm[i-1] ) && !isZero( relPerm[i] ) )
    {
      minVolFraction = phaseVolFrac[i-1];
    }
  }
  return minVolFraction;
}

TableRelativePermeability::KernelWrapper::
  KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & relPermKernelWrappers,
                 arrayView1d< real64 const > const & phaseMinVolumeFraction,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                 arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
  : RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac ),
  m_relPermKernelWrappers( relPermKernelWrappers ),
  m_phaseMinVolumeFraction( phaseMinVolumeFraction )
{}

TableRelativePermeability::KernelWrapper
TableRelativePermeability::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  // then we create the actual TableRelativePermeability::KernelWrapper
  return KernelWrapper( m_relPermKernelWrappers,
                        m_phaseMinVolumeFraction,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeability, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
