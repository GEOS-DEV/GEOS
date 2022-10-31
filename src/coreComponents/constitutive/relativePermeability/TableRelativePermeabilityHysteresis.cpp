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
 * @file TableRelativePermeabilityHysteresis.cpp
 */

#include "TableRelativePermeabilityHysteresis.hpp"

#include "constitutive/relativePermeability/RelativePermeabilityExtrinsicData.hpp"
#include "constitutive/relativePermeability/TableRelativePermeabilityHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableRelativePermeabilityHysteresis::TableRelativePermeabilityHysteresis( std::string const & name,
                                                                          Group * const parent )
  : RelativePermeabilityBase( name, parent )
{

  // drainage table names

  registerWrapper( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString(), &m_drainageWettingNonWettingRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of drainage relative permeability tables for the pair (wetting phase, non-wetting phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, nonWettingPhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString() ) +
                    " and " +
                    string( viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString(), &m_drainageWettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of drainage relative permeability tables for the pair (wetting phase, intermediate phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString(), &m_drainageNonWettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of drainage relative permeability tables for the pair (non-wetting phase, intermediate phase)\n"
                    "The expected format is \"{ nonWettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  // imbibition table names

  registerWrapper( viewKeyStruct::imbibitionWettingRelPermTableNameString(), &m_imbibitionWettingRelPermTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Imbibition relative permeability table name for the wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  registerWrapper( viewKeyStruct::imbibitionNonWettingRelPermTableNameString(), &m_imbibitionNonWettingRelPermTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Imbibition relative permeability table name for the non-wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  // hysteresis input parameters

  registerWrapper( viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "First parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Second parameter introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( viewKeyStruct::killoughCurvatureParameterString(), &m_killoughCurvatureParam ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Curvature parameter introduced by Killough for wetting-phase hysteresis (see RTD documentation)." );

  // internal class data

  registerWrapper( viewKeyStruct::drainagePhaseMinVolumeFractionString(), &m_drainagePhaseMinVolFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseMinVolumeFractionString(), &m_imbibitionPhaseMinVolFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::drainagePhaseRelPermEndPointString(), &m_drainagePhaseRelPermEndPoint ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseRelPermEndPointString(), &m_imbibitionPhaseRelPermEndPoint ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::drainagePhaseMaxVolumeFractionString(), &m_drainagePhaseMaxVolFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseMaxVolumeFractionString(), &m_imbibitionPhaseMaxVolFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::landParameterString(), &m_landParam ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseHasHysteresisString(), &m_phaseHasHysteresis ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerField( fields::relperm::phaseMaxHistoricalVolFraction{}, &m_phaseMaxHistoricalVolFraction );
  registerField( fields::relperm::phaseMinHistoricalVolFraction{}, &m_phaseMinHistoricalVolFraction );

  registerWrapper( viewKeyStruct::drainageRelPermKernelWrappersString(), &m_drainageRelPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::imbibitionRelPermKernelWrappersString(), &m_imbibitionRelPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );

}

void TableRelativePermeabilityHysteresis::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  integer const numPhases = m_phaseNames.size();
  GEOSX_THROW_IF( numPhases != 2 && numPhases != 3,
                  GEOSX_FMT( "{}: the expected number of fluid phases is either two, or three",
                             getFullName() ),
                  InputError );

  m_phaseHasHysteresis.resize( 2 );

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( m_drainageWettingNonWettingRelPermTableNames.empty(),
                    GEOSX_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the relative permeability tables "
                               "for the pair (wetting phase, non-wetting phase)",
                               getFullName(),
                               viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_drainageWettingNonWettingRelPermTableNames.size() != 2,
                    GEOSX_FMT( "{}: for a two-phase flow simulation, we must use {} to specify exactly two names: "
                               "first the name of the wetting phase relperm table, second the name on the non-wetting phase relperm table",
                               getFullName(),
                               viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ),
                    InputError );

    m_phaseHasHysteresis[IPT::WETTING] = ( m_imbibitionWettingRelPermTableName.empty() ||
                                           m_imbibitionWettingRelPermTableName == m_drainageWettingNonWettingRelPermTableNames[0] )
      ? 0 : 1;
    m_phaseHasHysteresis[IPT::NONWETTING] = ( m_imbibitionNonWettingRelPermTableName.empty() ||
                                              m_imbibitionNonWettingRelPermTableName == m_drainageWettingNonWettingRelPermTableNames[1] )
      ? 0 : 1;
  }
  else if( numPhases == 3 )
  {
    GEOSX_THROW_IF( m_drainageWettingIntermediateRelPermTableNames.empty() || m_drainageNonWettingIntermediateRelPermTableNames.empty(),
                    GEOSX_FMT( "{}: for a three-phase flow simulation, "
                               "we must use {} to specify the relative permeability tables for the pair (wetting phase, intermediate phase), "
                               "and {} to specify the relative permeability tables for the pair (non-wetting phase, intermediate phase)",
                               getFullName(),
                               viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString(),
                               viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString()  ),
                    InputError );

    GEOSX_THROW_IF( m_drainageWettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT( "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: "
                               "first the name of the wetting phase relperm table, second the name on the intermediate phase relperm table",
                               getFullName(),
                               viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_drainageNonWettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT( "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: "
                               "first the name of the non-wetting phase relperm table, second the name on the intermediate phase relperm table",
                               getFullName(),
                               viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString() ),
                    InputError );

    m_phaseHasHysteresis[IPT::WETTING] = ( m_imbibitionWettingRelPermTableName.empty() ||
                                           m_imbibitionWettingRelPermTableName == m_drainageWettingIntermediateRelPermTableNames[0] )
      ? 0 : 1;
    m_phaseHasHysteresis[IPT::NONWETTING] = ( m_imbibitionNonWettingRelPermTableName.empty() ||
                                              m_imbibitionNonWettingRelPermTableName == m_drainageNonWettingIntermediateRelPermTableNames[0] )
      ? 0 : 1;
  }

  GEOSX_THROW_IF( m_phaseHasHysteresis[IPT::WETTING] == 0 && m_phaseHasHysteresis[IPT::NONWETTING] == 0,
                  GEOSX_FMT( "{}: we must use {} or {} to specify at least one imbibition relative permeability table",
                             getFullName(),
                             viewKeyStruct::imbibitionWettingRelPermTableNameString(),
                             viewKeyStruct::imbibitionNonWettingRelPermTableNameString() ),
                  InputError );

  GEOSX_THROW_IF( m_jerauldParam_a < 0,
                  GEOSX_FMT( "{}: the parameter {} must be positive",
                             getFullName(),
                             viewKeyStruct::jerauldParameterAString() ),
                  InputError );

  GEOSX_THROW_IF( m_jerauldParam_b < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::jerauldParameterBString() ),
                  InputError );

  GEOSX_THROW_IF( m_killoughCurvatureParam < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::killoughCurvatureParameterString() ),
                  InputError );

}

void TableRelativePermeabilityHysteresis::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();

  m_drainagePhaseMinVolFraction.resize( MAX_NUM_PHASES );
  m_imbibitionPhaseMinVolFraction.resize( 2 ); // we don't save the value of the intermediate phase, for which we neglect hysteresis

  m_drainagePhaseMaxVolFraction.resize( MAX_NUM_PHASES );
  m_imbibitionPhaseMaxVolFraction.resize( 2 );

  m_drainagePhaseRelPermEndPoint.resize( MAX_NUM_PHASES );
  m_imbibitionPhaseRelPermEndPoint.resize( 2 ); // we don't save the value of the intermediate phase, for which we neglect hysteresis

  // Step 1: validate drainage relative permeabilities

  checkExistenceAndValidateDrainageRelPermTables();

  // Step 2: validate imbibition relative permeability tables

  checkExistenceAndValidateImbibitionRelPermTables();

  // Step 3: compute the Land coefficient

  computeLandCoefficient();
}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateDrainageRelPermTables()
{
  integer const numPhases = m_phaseNames.size();

  // Step 1.a: take care of the two-phase case

  if( numPhases == 2 )
  {
    for( integer ip = 0; ip < m_drainageWettingNonWettingRelPermTableNames.size(); ++ip )
    {
      if( ip == 0 ) // wetting phase is either water, or oil (for two-phase oil-gas systems)
      {
        integer const ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
        checkExistenceAndValidateRelPermTable( m_drainageWettingNonWettingRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[ipWetting], // output
                                               m_drainagePhaseMaxVolFraction[ipWetting],
                                               m_drainagePhaseRelPermEndPoint[ipWetting] );
      }
      else if( ip == 1 ) // non-wetting phase is either oil (for two-phase oil-water systems), or gas
      {
        integer const ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
        checkExistenceAndValidateRelPermTable( m_drainageWettingNonWettingRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[ipNonWetting], // output
                                               m_drainagePhaseMaxVolFraction[ipNonWetting],
                                               m_drainagePhaseRelPermEndPoint[ipNonWetting] );
      }
    }
  }
  // Step 1.b: take care of the three-phase case

  else if( numPhases == 3 )
  {
    for( integer ip = 0; ip < m_drainageWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      if( ip == 0 ) // wetting phase is water
      {
        checkExistenceAndValidateRelPermTable( m_drainageWettingIntermediateRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[m_phaseOrder[PhaseType::WATER]], // output
                                               m_drainagePhaseMaxVolFraction[m_phaseOrder[PhaseType::WATER]],
                                               m_drainagePhaseRelPermEndPoint[m_phaseOrder[PhaseType::WATER]] );
      }
      else if( ip == 1 ) // intermediate phase is oil
      {
        checkExistenceAndValidateRelPermTable( m_drainageWettingIntermediateRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[m_phaseOrder[PhaseType::OIL]], // output
                                               m_drainagePhaseMaxVolFraction[m_phaseOrder[PhaseType::OIL]],
                                               m_drainagePhaseRelPermEndPoint[m_phaseOrder[PhaseType::OIL]] );
      }
    }

    for( integer ip = 0; ip < m_drainageNonWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      if( ip == 0 ) // non-wetting phase is gas
      {
        checkExistenceAndValidateRelPermTable( m_drainageNonWettingIntermediateRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[m_phaseOrder[PhaseType::GAS]], // output
                                               m_drainagePhaseMaxVolFraction[m_phaseOrder[PhaseType::GAS]],
                                               m_drainagePhaseRelPermEndPoint[m_phaseOrder[PhaseType::GAS]] );
      }
      else if( ip == 1 ) // intermediate phase is oil
      {
        checkExistenceAndValidateRelPermTable( m_drainageNonWettingIntermediateRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[m_phaseOrder[PhaseType::OIL]], // output
                                               m_drainagePhaseMaxVolFraction[m_phaseOrder[PhaseType::OIL]],
                                               m_drainagePhaseRelPermEndPoint[m_phaseOrder[PhaseType::OIL]] );
      }
    }
  }
}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateImbibitionRelPermTables()
{
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  integer ipWetting = 0;
  integer ipNonWetting = 0;

  integer const numPhases = m_phaseNames.size();
  if( numPhases == 2 )
  {
    ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
    ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
  }
  else if( numPhases == 3 )
  {
    ipWetting = m_phaseOrder[PhaseType::WATER];
    ipNonWetting = m_phaseOrder[PhaseType::GAS];
  }

  // Step 1: validate wetting-phase imbibition relative permeability table

  if( m_phaseHasHysteresis[IPT::WETTING] )
  {

    checkExistenceAndValidateRelPermTable( m_imbibitionWettingRelPermTableName, // input
                                           m_imbibitionPhaseMinVolFraction[IPT::WETTING], // output
                                           m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                           m_imbibitionPhaseRelPermEndPoint[IPT::WETTING] );

    GEOSX_THROW_IF( !isZero( m_imbibitionPhaseMinVolFraction[IPT::WETTING] - m_drainagePhaseMinVolFraction[ipWetting] ),
                    GEOSX_FMT( "{}: the critical wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n"
                               "However, we found that the drainage critical wetting-phase volume fraction is {}, "
                               "whereas the imbibition critical wetting-phase volume fraction is {}",
                               getFullName(),
                               m_drainagePhaseMinVolFraction[ipWetting], m_imbibitionPhaseMinVolFraction[IPT::WETTING] ),
                    InputError );

    GEOSX_THROW_IF( m_imbibitionPhaseMaxVolFraction[IPT::WETTING] > m_drainagePhaseMaxVolFraction[ipWetting],
                    GEOSX_FMT( "{}: the maximum wetting-phase volume fraction (saturation) must be smaller in imbibition (compared to the drainage value).\n"
                               "However, we found that the drainage maximum wetting-phase volume fraction is {}, "
                               "whereas the imbibition maximum wetting-phase volume fraction is {}",
                               getFullName(),
                               m_drainagePhaseMaxVolFraction[ipWetting], m_imbibitionPhaseMaxVolFraction[IPT::WETTING] ),
                    InputError );
  }

  // Step 2: validate non-wetting-phase imbibition relative permeability table

  if( m_phaseHasHysteresis[IPT::NONWETTING] )
  {

    checkExistenceAndValidateRelPermTable( m_imbibitionNonWettingRelPermTableName, // input
                                           m_imbibitionPhaseMinVolFraction[IPT::NONWETTING], // output
                                           m_imbibitionPhaseMaxVolFraction[IPT::NONWETTING],
                                           m_imbibitionPhaseRelPermEndPoint[IPT::NONWETTING] );

    GEOSX_THROW_IF( !isZero ( m_imbibitionPhaseMaxVolFraction[IPT::NONWETTING] - m_drainagePhaseMaxVolFraction[ipNonWetting] ),
                    GEOSX_FMT( string( "{}: the maximum non-wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage maximum wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition maximum wetting-phase volume fraction is {}" ),
                               getFullName(),
                               m_drainagePhaseMaxVolFraction[ipNonWetting], m_imbibitionPhaseMaxVolFraction[IPT::NONWETTING] ),
                    InputError );

    GEOSX_THROW_IF( !isZero ( m_imbibitionPhaseRelPermEndPoint[IPT::NONWETTING] - m_drainagePhaseRelPermEndPoint[ipNonWetting] ),
                    GEOSX_FMT( string( "{}: the non-wetting-phase relperm endpoint must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage endpoint wetting-phase relperm is {}, " )
                               + string( "whereas the imbibition endpoint wetting-phase relperm is {}" ),
                               getFullName(),
                               m_drainagePhaseRelPermEndPoint[ipNonWetting], m_imbibitionPhaseRelPermEndPoint[IPT::NONWETTING] ),
                    InputError );

    GEOSX_THROW_IF( m_imbibitionPhaseMinVolFraction[IPT::NONWETTING] < m_drainagePhaseMinVolFraction[ipNonWetting],
                    GEOSX_FMT( string( "{}: the critical wetting-phase volume fraction (saturation) must be larger in imbibition (compared to the drainage value).\n" )
                               + string( "However, we found that the drainage critical wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition critical wetting-phase volume fraction is {}" ),
                               getFullName(),
                               m_drainagePhaseMinVolFraction[ipNonWetting], m_imbibitionPhaseMinVolFraction[IPT::NONWETTING] ),
                    InputError );
  }
}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateRelPermTable( string const & relPermTableName,
                                                                                 real64 & phaseMinVolFrac,
                                                                                 real64 & phaseMaxVolFrac,
                                                                                 real64 & phaseRelPermEndPoint ) const
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  // check if the table actually exists
  GEOSX_THROW_IF( !functionManager.hasGroup( relPermTableName ),
                  GEOSX_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             relPermTableName ),
                  InputError );
  TableFunction const & relPermTable = functionManager.getGroup< TableFunction >( relPermTableName );

  // read the table, check monotonicity, and return the min/max saturation and the endpoint
  string const fullName = getFullName();
  TableRelativePermeabilityHelpers::validateRelativePermeabilityTable( relPermTable, // input
                                                                       fullName,
                                                                       phaseMinVolFrac, // output
                                                                       phaseMaxVolFrac,
                                                                       phaseRelPermEndPoint );
}

void TableRelativePermeabilityHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

  integer ipWetting = 0;
  integer ipNonWetting = 0;

  integer const numPhases = m_phaseNames.size();
  if( numPhases == 2 )
  {
    ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
    ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
  }
  else if( numPhases == 3 )
  {
    ipWetting = m_phaseOrder[PhaseType::WATER];
    ipNonWetting = m_phaseOrder[PhaseType::GAS];
  }

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase

  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  {
    real64 const Scrd = m_drainagePhaseMinVolFraction[ipWetting];
    real64 const Smxd = m_drainagePhaseMaxVolFraction[ipWetting];
    real64 const Smxi = m_imbibitionPhaseMaxVolFraction[IPT::WETTING];
    real64 const Swc = Scrd;
    GEOSX_THROW_IF(  (Smxi - Smxd) > 0,
                     GEOSX_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                                "Crossing relative permeability curves.\n",
                                getFullName(),
                                Smxi,
                                Smxd ),
                     InputError );

    m_landParam[IPT::WETTING] = ( Smxd - Swc ) / LvArray::math::max( KernelWrapper::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
  }

  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Scrd = m_drainagePhaseMinVolFraction[ipNonWetting];
    real64 const Scri = m_imbibitionPhaseMinVolFraction[IPT::NONWETTING];
    real64 const Smx = m_drainagePhaseMaxVolFraction[ipNonWetting];
    GEOSX_THROW_IF( (Scrd - Scri) > 0,
                    GEOSX_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                               "Crossing relative permeability curves.\n",
                               getFullName(),
                               Scrd,
                               Scri ),
                    InputError );

    m_landParam[IPT::NONWETTING] = ( Smx - Scrd ) / LvArray::math::max( KernelWrapper::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
  }
}

void TableRelativePermeabilityHysteresis::createAllTableKernelWrappers()
{
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_drainageRelPermKernelWrappers.clear();
  m_imbibitionRelPermKernelWrappers.clear();

  if( numPhases == 2 )
  {
    for( integer ip = 0; ip < m_drainageWettingNonWettingRelPermTableNames.size(); ++ip )
    {
      TableFunction const & drainageRelPermTable = functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingRelPermTableNames[ip] );
      m_drainageRelPermKernelWrappers.emplace_back( drainageRelPermTable.createKernelWrapper() );
    }

    TableFunction const & imbibitionWettingRelPermTable = m_phaseHasHysteresis[IPT::WETTING]
      ? functionManager.getGroup< TableFunction >( m_imbibitionWettingRelPermTableName )
      : functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingRelPermTableNames[0] );
    m_imbibitionRelPermKernelWrappers.emplace_back( imbibitionWettingRelPermTable.createKernelWrapper() );

    TableFunction const & imbibitionNonWettingRelPermTable = m_phaseHasHysteresis[IPT::NONWETTING]
      ? functionManager.getGroup< TableFunction >( m_imbibitionNonWettingRelPermTableName )
      : functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingRelPermTableNames[1] );
    m_imbibitionRelPermKernelWrappers.emplace_back( imbibitionNonWettingRelPermTable.createKernelWrapper() );

  }
  else if( numPhases == 3 )
  {
    for( integer ip = 0; ip < m_drainageWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      TableFunction const & drainageRelPermTable = functionManager.getGroup< TableFunction >( m_drainageWettingIntermediateRelPermTableNames[ip] );
      m_drainageRelPermKernelWrappers.emplace_back( drainageRelPermTable.createKernelWrapper() );
    }
    for( integer ip = 0; ip < m_drainageNonWettingIntermediateRelPermTableNames.size(); ++ip )
    {
      TableFunction const & drainageRelPermTable = functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateRelPermTableNames[ip] );
      m_drainageRelPermKernelWrappers.emplace_back( drainageRelPermTable.createKernelWrapper() );
    }

    TableFunction const & imbibitionWettingRelPermTable = m_phaseHasHysteresis[IPT::WETTING]
      ? functionManager.getGroup< TableFunction >( m_imbibitionWettingRelPermTableName )
      : functionManager.getGroup< TableFunction >( m_drainageWettingIntermediateRelPermTableNames[0] );
    m_imbibitionRelPermKernelWrappers.emplace_back( imbibitionWettingRelPermTable.createKernelWrapper() );

    TableFunction const & imbibitionNonWettingRelPermTable = m_phaseHasHysteresis[IPT::NONWETTING]
      ? functionManager.getGroup< TableFunction >( m_imbibitionNonWettingRelPermTableName )
      : functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateRelPermTableNames[0] );
    m_imbibitionRelPermKernelWrappers.emplace_back( imbibitionNonWettingRelPermTable.createKernelWrapper() );
  }

}

TableRelativePermeabilityHysteresis::KernelWrapper
TableRelativePermeabilityHysteresis::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
  return KernelWrapper( m_drainageRelPermKernelWrappers,
                        m_imbibitionRelPermKernelWrappers,
                        m_jerauldParam_a,
                        m_jerauldParam_b,
                        m_killoughCurvatureParam,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_drainagePhaseMinVolFraction,
                        m_imbibitionPhaseMinVolFraction,
                        m_drainagePhaseMaxVolFraction,
                        m_imbibitionPhaseMaxVolFraction,
                        m_drainagePhaseRelPermEndPoint,
                        m_imbibitionPhaseRelPermEndPoint,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseMinHistoricalVolFraction,
                        m_phaseMaxHistoricalVolFraction,
                        m_phaseTrappedVolFrac,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

void TableRelativePermeabilityHysteresis::resizeFields( localIndex const size, localIndex const numPts )
{
  RelativePermeabilityBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();

  m_phaseMaxHistoricalVolFraction.resize( size, numPhases );
  m_phaseMinHistoricalVolFraction.resize( size, numPhases );
  m_phaseMaxHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
  m_phaseMinHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 1.0 );
}

void TableRelativePermeabilityHysteresis::saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
{
  RelativePermeabilityBase::saveConvergedState();

  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();

  localIndex const numElems = phaseVolFraction.size( 0 );
  integer const numPhases = numFluidPhases();

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMaxHistoricalVolFraction[ei][ip] = LvArray::math::max( phaseVolFraction[ei][ip], phaseMaxHistoricalVolFraction[ei][ip] );
      phaseMinHistoricalVolFraction[ei][ip] = LvArray::math::min( phaseVolFraction[ei][ip], phaseMinHistoricalVolFraction[ei][ip] );
    }
  } );

}

TableRelativePermeabilityHysteresis::KernelWrapper::KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & drainageRelPermKernelWrappers,
                                                                   arrayView1d< TableFunction::KernelWrapper const > const & imbibitionRelPermKernelWrappers,
                                                                   real64 const & jerauldParam_a,
                                                                   real64 const & jerauldParam_b,
                                                                   real64 const & killoughCurvatureParam,
                                                                   arrayView1d< integer const > const & phaseHasHysteresis,
                                                                   arrayView1d< real64 const > const & landParam,
                                                                   arrayView1d< real64 const > const & drainagePhaseMinVolFraction,
                                                                   arrayView1d< real64 const > const & imbibitionPhaseMinVolFraction,
                                                                   arrayView1d< real64 const > const & drainagePhaseMaxVolFraction,
                                                                   arrayView1d< real64 const > const & imbibitionPhaseMaxVolFraction,
                                                                   arrayView1d< real64 const > const & drainagePhaseRelPermEndPoint,
                                                                   arrayView1d< real64 const > const & imbibitionPhaseRelPermEndPoint,
                                                                   arrayView1d< integer const > const & phaseTypes,
                                                                   arrayView1d< integer const > const & phaseOrder,
                                                                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
                                                                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
                                                                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrappedVolFrac,
                                                                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                                                   arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
  : RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac,
                                    phaseTrappedVolFrac ),
  m_drainageRelPermKernelWrappers( drainageRelPermKernelWrappers ),
  m_imbibitionRelPermKernelWrappers( imbibitionRelPermKernelWrappers ),
  m_jerauldParam_a( jerauldParam_a ),
  m_jerauldParam_b( jerauldParam_b ),
  m_killoughCurvatureParam( killoughCurvatureParam ),
  m_phaseHasHysteresis( phaseHasHysteresis ),
  m_landParam( landParam ),
  m_drainagePhaseMinVolFraction( drainagePhaseMinVolFraction ),
  m_imbibitionPhaseMinVolFraction( imbibitionPhaseMinVolFraction ),
  m_drainagePhaseMaxVolFraction( drainagePhaseMaxVolFraction ),
  m_imbibitionPhaseMaxVolFraction( imbibitionPhaseMaxVolFraction ),
  m_drainagePhaseRelPermEndPoint( drainagePhaseRelPermEndPoint ),
  m_imbibitionPhaseRelPermEndPoint( imbibitionPhaseRelPermEndPoint ),
  m_phaseMinHistoricalVolFraction( phaseMinHistoricalVolFraction ),
  m_phaseMaxHistoricalVolFraction( phaseMaxHistoricalVolFraction )
{}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeabilityHysteresis, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
