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
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Imbibition relative permeability table name for the wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  registerWrapper( viewKeyStruct::imbibitionNonWettingRelPermTableNameString(), &m_imbibitionNonWettingRelPermTableName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Imbibition relative permeability table name for the non-wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  // hysteresis input parameters

  registerWrapper( viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "First parameter introduced by Jerauld (see RTD documentation)." );

  registerWrapper( viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Second parameter introduced by Jerauld (see RTD documentation)." );

  registerWrapper( viewKeyStruct::alphaParameter2String(), &m_alphaParam_2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Second parameter introduced by Jerauld (see RTD documentation)." );

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

  localIndex const numPhases = m_phaseNames.size();
  GEOSX_THROW_IF( numPhases != 2 && numPhases != 3,
                  GEOSX_FMT( "{}: the expected number of fluid phases is either two, or three",
                             getFullName() ),
                  InputError );

  m_phaseHasHysteresis.resize( 2 );

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( m_drainageWettingNonWettingRelPermTableNames.empty(),
                    GEOSX_FMT( string( "{}: for a two-phase flow simulation, we must use {} to specify the relative permeability tables " )
                               + string( "for the pair (wetting phase, non-wetting phase)" ),
                               getFullName(),
                               viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_drainageWettingNonWettingRelPermTableNames.size() != 2,
                    GEOSX_FMT( string( "{}: for a two-phase flow simulation, we must use {} to specify exactly two names: " )
                               + string( "first the name of the wetting phase relperm table, second the name on the non-wetting phase relperm table" ),
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
                    GEOSX_FMT( string( "{}: for a three-phase flow simulation, " )
                               + string( "we must use {} to specify the relative permeability tables for the pair (wetting phase, intermediate phase), " )
                               + string( "and {} to specify the relative permeability tables for the pair (non-wetting phase, intermediate phase)" ),
                               getFullName(),
                               viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString(),
                               viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString()  ),
                    InputError );

    GEOSX_THROW_IF( m_drainageWettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT( string( "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: " )
                               + string( "first the name of the wetting phase relperm table, second the name on the intermediate phase relperm table" ),
                               getFullName(),
                               viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString() ),
                    InputError );

    GEOSX_THROW_IF( m_drainageNonWettingIntermediateRelPermTableNames.size() != 2,
                    GEOSX_FMT( string( "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: " )
                               + string( "first the name of the non-wetting phase relperm table, second the name on the intermediate phase relperm table" ),
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

  GEOSX_THROW_IF( m_alphaParam_2 < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::alphaParameter2String() ),
                  InputError );

}

void TableRelativePermeabilityHysteresis::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();

  m_drainagePhaseMinVolFraction.resize( MAX_NUM_PHASES );
  m_imbibitionPhaseMinVolFraction.resize( 2 ); // we don't save the value of the intermediate phase, for which we neglect hysteresis

  m_drainagePhaseMaxVolFraction.resize( MAX_NUM_PHASES );// two values as we consider case of wetting phase hysteresis
  m_imbibitionPhaseMaxVolFraction.resize( MAX_NUM_PHASES );

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
  localIndex const numPhases = m_phaseNames.size();

  // Step 1.a: take care of the two-phase case

  if( numPhases == 2 )
  {
    for( integer ip = 0; ip < m_drainageWettingNonWettingRelPermTableNames.size(); ++ip )
    {
      if( ip == 0 ) // wetting phase is either water, or oil (for two-phase oil-gas systems)
      {
        localIndex const ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
        checkExistenceAndValidateRelPermTable( m_drainageWettingNonWettingRelPermTableNames[ip], // input
                                               m_drainagePhaseMinVolFraction[ipWetting], // output
                                               m_drainagePhaseMaxVolFraction[ipWetting],
                                               m_drainagePhaseRelPermEndPoint[ipWetting] );
      }
      else if( ip == 1 ) // non-wetting phase is either oil (for two-phase oil-water systems), or gas
      {
        localIndex const ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
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

  localIndex ipWetting = 0;
  localIndex ipNonWetting = 0;
//  real64 maxVolFraction = 0;

  localIndex const numPhases = m_phaseNames.size();
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

  // Step 2.a: validate wetting-phase imbibition relative permeability table

  if( m_phaseHasHysteresis[IPT::WETTING] )
  {

    checkExistenceAndValidateRelPermTable( m_imbibitionWettingRelPermTableName, // input
                                           m_imbibitionPhaseMinVolFraction[IPT::WETTING], // output
                                           m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                           m_imbibitionPhaseRelPermEndPoint[IPT::WETTING] );

    GEOSX_THROW_IF( !isZero( m_imbibitionPhaseMinVolFraction[IPT::WETTING] - m_drainagePhaseMinVolFraction[ipWetting] ),
                    GEOSX_FMT( string( "{}: the critical wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage critical wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition critical wetting-phase volume fraction is {}" ),
                               getFullName(),
                               m_drainagePhaseMinVolFraction[ipWetting], m_imbibitionPhaseMinVolFraction[IPT::WETTING] ),
                    InputError );

    GEOSX_THROW_IF( m_imbibitionPhaseMaxVolFraction[ipWetting] > m_drainagePhaseMaxVolFraction[ipWetting],
                    GEOSX_FMT( string( "{}: the maximum wetting-phase volume fraction (saturation) must be smaller in imbibition (compared to the drainage value).\n" )
                               + string( "However, we found that the drainage maximum wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition maximum wetting-phase volume fraction is {}" ),
                               getFullName(),
                               m_drainagePhaseMaxVolFraction[ipWetting], m_imbibitionPhaseMaxVolFraction[ipWetting]),
                    InputError );
  }

  // Step 2.b: validate non-wetting-phase imbibition relative permeability table

  if( m_phaseHasHysteresis[IPT::NONWETTING] )
  {

    checkExistenceAndValidateRelPermTable( m_imbibitionNonWettingRelPermTableName, // input
                                           m_imbibitionPhaseMinVolFraction[IPT::NONWETTING], // output
                                           m_imbibitionPhaseMaxVolFraction[ipNonWetting],
                                           m_imbibitionPhaseRelPermEndPoint[IPT::NONWETTING] );

    GEOSX_THROW_IF( !isZero ( m_imbibitionPhaseMaxVolFraction[ipNonWetting] - m_drainagePhaseMaxVolFraction[ipNonWetting] ),
                    GEOSX_FMT( string( "{}: the maximum non-wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage maximum wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition maximum wetting-phase volume fraction is {}" ),
                               getFullName(),
                               m_drainagePhaseMaxVolFraction[ipNonWetting] , m_imbibitionPhaseMaxVolFraction[ipNonWetting] ),
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
  validateRelPermTable( relPermTable, // input
                        phaseMinVolFrac, // output
                        phaseMaxVolFrac,
                        phaseRelPermEndPoint );
}

void TableRelativePermeabilityHysteresis::validateRelPermTable( TableFunction const & relPermTable,
                                                                real64 & phaseMinVolFrac,
                                                                real64 & phaseMaxVolFrac,
                                                                real64 & phaseRelPermEndPoint ) const
{
  // TODO: merge this function with TableRelativePermeability::validateRelPermTable

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
  phaseMinVolFrac = phaseVolFrac[0];
  phaseMaxVolFrac = phaseVolFrac[phaseVolFrac.size()-1];
  phaseRelPermEndPoint = relPerm[relPerm.size()-1];

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
//    GEOSX_THROW_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-10,
    GEOSX_THROW_IF(  (relPerm[i] - relPerm[i-1]) < -1e-10,
                    GEOSX_FMT( "{}: in table '{}' values must be strictly increasing", getFullName(), relPermTable.getName() ),
                    InputError );

    if( isZero( relPerm[i-1] ) && !isZero( relPerm[i] ) )
    {
      phaseMinVolFrac = phaseVolFrac[i-1];
    }
    if( !isZero( relPerm[i-1] ) && isZero( relPerm[i] - relPerm[i-1] ) )
    {
      phaseMaxVolFrac = phaseVolFrac[i-1];
    }
  }

}

void TableRelativePermeabilityHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

  localIndex ipWetting = 0;
  localIndex ipNonWetting = 0;

  localIndex const numPhases = m_phaseNames.size();
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

    real64 const Smxd = m_drainagePhaseMaxVolFraction[ipWetting];
    real64 const Smxi = m_imbibitionPhaseMaxVolFraction[ipWetting];
    real64 const Scrd = m_drainagePhaseMinVolFraction[IPT::WETTING];
    real64 const Scri = m_imbibitionPhaseMinVolFraction[IPT::WETTING];
    real64 const Swc = Scrd;
    GEOSX_ERROR_IF(  (Smxi - Smxd) > 0,
     GEOSX_FMT( string("{}: For wetting phase hysteresis, imbibition end-point Smxi( {} ) is larger than drainage one Smxd( {} ).\n") +
                string("Crossing curves.\n"),
                getFullName(),
                Smxi,
                Smxd )
      );

    m_landParam[IPT::WETTING] = ( Smxd - Swc ) / ( Smxd - Smxi ) - 1.0;
  }

  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Smx = m_drainagePhaseMaxVolFraction[ipNonWetting];
    real64 const Scrd = m_drainagePhaseMinVolFraction[ipNonWetting];
    real64 const Scri = m_imbibitionPhaseMinVolFraction[IPT::NONWETTING];
    GEOSX_ERROR_IF( (Scrd - Scri) > 0,
      GEOSX_FMT(string("{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) is larger than imbibition one Scri( {} ).\n") +
                string("Crossing curves.\n"),
                getFullName(),
                Scrd,
                Scri )
    );

    m_landParam[IPT::NONWETTING] = ( Smx - Scrd ) / ( Scri - Scrd ) - 1.0;
  }

  // Step 3: make sure that they match for two-phase flow

  if( m_phaseHasHysteresis[IPT::WETTING] && m_phaseHasHysteresis[IPT::NONWETTING] )
  {
    GEOSX_WARNING_IF( numPhases == 2 && !isZero( m_landParam[IPT::WETTING] - m_landParam[IPT::NONWETTING] ),
                    GEOSX_FMT( string( "{}: For two-phase flow, the Land parameters computed from the wetting and non-wetting relperm curves must match.\n" )
                               + string( "However, we found that the wetting Land parameter is {}, " )
                               + string( "whereas the nonwetting Land parameter is {}. " )
                               + string( "This might result in inconsistency." ),
                               getFullName(),
                               m_landParam[IPT::WETTING], m_landParam[IPT::NONWETTING] ) );
  }
}

void TableRelativePermeabilityHysteresis::createAllTableKernelWrappers()
{
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  FunctionManager const & functionManager = FunctionManager::getInstance();

  localIndex const numPhases = m_phaseNames.size();

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
                        m_alphaParam_2,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_drainagePhaseMinVolFraction,
                        m_imbibitionPhaseMinVolFraction,
                        m_drainagePhaseMaxVolFraction,
                        m_imbibitionPhaseMaxVolFraction,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseMaxHistoricalVolFraction,
                        m_phaseMinHistoricalVolFraction,
                        m_drainagePhaseRelPermEndPoint,
                        m_imbibitionPhaseRelPermEndPoint,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

void TableRelativePermeabilityHysteresis::resizeFields( localIndex const size, localIndex const numPts )
{
  RelativePermeabilityBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();

  m_phaseMaxHistoricalVolFraction.resize( size, numPhases );
  m_phaseMinHistoricalVolFraction.resize( size, numPhases );
}

void TableRelativePermeabilityHysteresis::initializeState( arrayView2d< real64 const, compflow::USD_PHASE > const & initialPhaseVolFraction ) const
{
  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();

  localIndex const numElems = initialPhaseVolFraction.size( 0 );
  localIndex const numPhases = numFluidPhases();

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMaxHistoricalVolFraction[ei][ip] = initialPhaseVolFraction[ei][ip];
      phaseMinHistoricalVolFraction[ei][ip] = initialPhaseVolFraction[ei][ip];
    }
  } );
}

void TableRelativePermeabilityHysteresis::saveConvergedPhaseVolFraction( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
{
  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();

  localIndex const numElems = phaseVolFraction.size( 0 );
  localIndex const numPhases = numFluidPhases();

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      if( phaseVolFraction[ei][ip] > phaseMaxHistoricalVolFraction[ei][ip] )
      {
        phaseMaxHistoricalVolFraction[ei][ip] = phaseVolFraction[ei][ip];
      }
      if( phaseVolFraction[ei][ip] < phaseMinHistoricalVolFraction[ei][ip] )
      {
        phaseMinHistoricalVolFraction[ei][ip] = phaseVolFraction[ei][ip];
      }
    }
  } );

}

TableRelativePermeabilityHysteresis::KernelWrapper::
  KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & drainageRelPermKernelWrappers,
                 arrayView1d< TableFunction::KernelWrapper const > const & imbibitionRelPermKernelWrappers,
                 real64 const & jerauldParam_a,
                 real64 const & jerauldParam_b,
                 real64 const & alphaParam_2,
                 arrayView1d< integer const > const & phaseHasHysteresis,
                 arrayView1d< real64 const > const & landParam,
                 arrayView1d< real64 const > const & drainagePhaseMinVolFraction,
                 arrayView1d< real64 const > const & imbibitionPhaseMinVolFraction,
                 arrayView1d< real64 const > const & drainagePhaseMaxVolFraction,
                 arrayView1d< real64 const > const & imbibitionPhaseMaxVolFraction,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
                 arrayView1d< real64 const > const & drainagePhaseRelPermEndPoint,
                 arrayView1d< real64 const > const & imbibitionPhaseRelPermEndPoint,
                 arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                 arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
  : RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac ),
  m_drainageRelPermKernelWrappers( drainageRelPermKernelWrappers ),
  m_imbibitionRelPermKernelWrappers( imbibitionRelPermKernelWrappers ),
  m_jerauldParam_a( jerauldParam_a ),
  m_jerauldParam_b( jerauldParam_b ),
  m_alphaParam_2( alphaParam_2 ),
  m_phaseHasHysteresis( phaseHasHysteresis ),
  m_landParam( landParam ),
  m_drainagePhaseMinVolFraction( drainagePhaseMinVolFraction ),
  m_imbibitionPhaseMinVolFraction( imbibitionPhaseMinVolFraction ),
  m_drainagePhaseMaxVolFraction( drainagePhaseMaxVolFraction ),
  m_imbibitionPhaseMaxVolFraction( imbibitionPhaseMaxVolFraction ),
  m_phaseMaxHistoricalVolFraction( phaseMaxHistoricalVolFraction ),
  m_phaseMinHistoricalVolFraction( phaseMinHistoricalVolFraction ),
  m_drainagePhaseRelPermEndPoint( drainagePhaseRelPermEndPoint ),
  m_imbibitionPhaseRelPermEndPoint( imbibitionPhaseRelPermEndPoint )
{}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeabilityHysteresis, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
