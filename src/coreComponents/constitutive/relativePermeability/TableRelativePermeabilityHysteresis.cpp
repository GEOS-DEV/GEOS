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

#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
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

  registerWrapper( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString(),
                   &m_drainageWettingNonWettingRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "List of drainage relative permeability tables for the pair (wetting phase, non-wetting phase)\n"
    "The expected format is \"{ wettingPhaseRelPermTableName, nonWettingPhaseRelPermTableName }\", in that order\n"
    "Note that this input is only used for two-phase flow.\n"
    "If you want to do a three-phase simulation, please use instead " +
    string( viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString() ) +
    " and " +
    string( viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString() ) +
    " to specify the table names" );

  registerWrapper( viewKeyStruct::drainageWettingIntermediateRelPermTableNamesString(),
                   &m_drainageWettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "List of drainage relative permeability tables for the pair (wetting phase, intermediate phase)\n"
    "The expected format is \"{ wettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
    "Note that this input is only used for three-phase flow.\n"
    "If you want to do a two-phase simulation, please use instead " +
    string( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ) +
    " to specify the table names" );

  registerWrapper( viewKeyStruct::drainageNonWettingIntermediateRelPermTableNamesString(),
                   &m_drainageNonWettingIntermediateRelPermTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "List of drainage relative permeability tables for the pair (non-wetting phase, intermediate phase)\n"
    "The expected format is \"{ nonWettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
    "Note that this input is only used for three-phase flow.\n"
    "If you want to do a two-phase simulation, please use instead " +
    string( viewKeyStruct::drainageWettingNonWettingRelPermTableNamesString() ) +
    " to specify the table names" );

  // imbibition table names

  registerWrapper( viewKeyStruct::imbibitionWettingRelPermTableNameString(), &m_imbibitionWettingRelPermTableName ).
    setInputFlag(
    InputFlags::OPTIONAL ).
    setApplyDefaultValue(
    "" ).
    setDescription(
    "Imbibition relative permeability table name for the wetting phase.\n"
    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  registerWrapper( viewKeyStruct::imbibitionNonWettingRelPermTableNameString(),
                   &m_imbibitionNonWettingRelPermTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription(
    "Imbibition relative permeability table name for the non-wetting phase.\n"
    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  // hysteresis input parameters
  registerWrapper( viewKeyStruct::phaseHasHysteresisString(), &m_phaseHasHysteresis ).
    setInputFlag( InputFlags::FALSE )
    .                                                                                  // will be deduced from tables
    setSizedFromParent( 0 );

  registerField( fields::relperm::phaseMaxHistoricalVolFraction{},
                 &m_phaseMaxHistoricalVolFraction );
  registerField( fields::relperm::phaseMinHistoricalVolFraction{},
                 &m_phaseMinHistoricalVolFraction );

  /// Killough data
  registerWrapper( viewKeyStruct::KilloughModelString(), &m_KilloughModel )
    .setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );


  registerWrapper( viewKeyStruct::KilloughModelWrapperString(), &m_KilloughKernel )
    .setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );


  registerWrapper( viewKeyStruct::drainageRelPermKernelWrappersString(),
                   &m_drainageRelPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags(
    RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::imbibitionRelPermKernelWrappersString(),
                   &m_imbibitionRelPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags(
    RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::landParameterString(), &m_landParam ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );


  registerWrapper( viewKeyStruct::wettingCurveString(), &m_wettingCurve ).
    setInputFlag(
    InputFlags::FALSE ).         // will be deduced from tables
    setSizedFromParent(
    0 )
    .setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::nonWettingCurveString(), &m_nonWettingCurve ).
    setInputFlag(
    InputFlags::FALSE ).         // will be deduced from tables
    setSizedFromParent(
    0 )
    .setRestartFlags( RestartFlags::NO_WRITE );

  //Forwarded to KilloughHysteresis
  registerWrapper( KilloughHysteresis::viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription(
    "First parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( KilloughHysteresis::viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription(
    "Second parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( KilloughHysteresis::viewKeyStruct::killoughCurvatureParameterString(), &m_killoughCurvatureParamRelPerm ).
    setInputFlag(
    InputFlags::OPTIONAL ).
    setApplyDefaultValue(
    1.0 ).
    setDescription(
    "Curvature parameter introduced by Killough for wetting-phase hysteresis (see RTD documentation)." );

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

  //Killough section
  m_KilloughModel.setRelPermParameters( m_jerauldParam_a, m_jerauldParam_b, m_killoughCurvatureParamRelPerm );
  m_KilloughModel.postProcessInput();
}

void TableRelativePermeabilityHysteresis::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();
  // Step 1: validate drainage relative permeabilities

  checkExistenceAndValidateWettingRelPermTables();

  // Step 2: validate imbibition relative permeability tables

  checkExistenceAndValidateNonWettingRelPermTables();

  //Step 3: validate intermediate permeability tables

  checkExistenceAndValidateIntermediateRelPermTables();

  // Step 4: compute the Land coefficient
  computeLandCoefficient();
}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateWettingRelPermTables()
{
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;
  integer const numPhases = m_phaseNames.size();
  integer ipWetting = -1, ipNonWetting = -1;
  std::tie( ipWetting, ipNonWetting ) = RelativePermeabilityBase::phaseIndex( m_phaseOrder );

  // Step 1.a: take care of the two-phase case
  real64 drainagePhaseMinVolFraction,     // output
         drainagePhaseMaxVolFraction,
         drainagePhaseRelPermMinEndPoint,
         drainagePhaseRelPermMaxEndPoint;
  GEOSX_ASSERT( m_drainageWettingNonWettingRelPermTableNames.size() == 2 );
  auto tableName = ( numPhases == 2 ) ?   m_drainageWettingNonWettingRelPermTableNames[0] : m_drainageWettingIntermediateRelPermTableNames[0];
//        integer const ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
  checkExistenceAndValidateRelPermTable( tableName,    // input
                                         drainagePhaseMinVolFraction,    // output
                                         drainagePhaseMaxVolFraction,
                                         drainagePhaseRelPermMinEndPoint,
                                         drainagePhaseRelPermMaxEndPoint );


  //imbibition if provided
  real64 imbibitionPhaseMinVolFraction,    // output
         imbibitionPhaseMaxVolFraction,
         imbibitionPhaseRelPermMinEndPoint,
         imbibitionPhaseRelPermMaxEndPoint;

  if( m_phaseHasHysteresis[IPT::WETTING] )
  {
    checkExistenceAndValidateRelPermTable( m_imbibitionWettingRelPermTableName,    // input
                                           imbibitionPhaseMinVolFraction,    // output
                                           imbibitionPhaseMaxVolFraction,
                                           imbibitionPhaseRelPermMinEndPoint,
                                           imbibitionPhaseRelPermMaxEndPoint );

    GEOSX_THROW_IF( !isZero( imbibitionPhaseMinVolFraction - drainagePhaseMinVolFraction ),
                    GEOSX_FMT( "{}: the critical wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n"
                               "However, we found that the drainage critical wetting-phase volume fraction is {}, "
                               "whereas the imbibition critical wetting-phase volume fraction is {}",
                               getFullName(),
                               drainagePhaseMinVolFraction, imbibitionPhaseMinVolFraction ),
                    InputError );

    GEOSX_THROW_IF( imbibitionPhaseMaxVolFraction > drainagePhaseMaxVolFraction,
                    GEOSX_FMT( "{}: the maximum wetting-phase volume fraction (saturation) must be smaller in imbibition (compared to the drainage value).\n"
                               "However, we found that the drainage maximum wetting-phase volume fraction is {}, "
                               "whereas the imbibition maximum wetting-phase volume fraction is {}",
                               getFullName(),
                               drainagePhaseMaxVolFraction, imbibitionPhaseMaxVolFraction ),
                    InputError );

    GEOSX_THROW_IF( imbibitionPhaseRelPermMaxEndPoint > drainagePhaseRelPermMaxEndPoint,
                    GEOSX_FMT( "{}: the maximum wetting-phase relperm must be smaller in imbibition (compared to the drainage value).\n"
                               "However, we found that the drainage maximum wetting-phase relperm is {}, "
                               "whereas the imbibition maximum wetting-phase relperm is {}",
                               getFullName(),
                               drainagePhaseRelPermMaxEndPoint, imbibitionPhaseRelPermMaxEndPoint ),
                    InputError );

    //test other end point ??

  }

  m_wettingCurve.setPoints( {drainagePhaseMinVolFraction, drainagePhaseRelPermMinEndPoint},   // same as imbibition min
                            {imbibitionPhaseMaxVolFraction, imbibitionPhaseRelPermMaxEndPoint},
                            {drainagePhaseMaxVolFraction, drainagePhaseRelPermMaxEndPoint} );
}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateNonWettingRelPermTables()
{
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  integer const numPhases = m_phaseNames.size();
  integer ipWetting = -1, ipNonWetting = -1;
  std::tie( ipWetting, ipNonWetting ) = RelativePermeabilityBase::phaseIndex( m_phaseOrder );

  //treat drainage
  real64 drainagePhaseMinVolFraction,   // output
         drainagePhaseMaxVolFraction,
         drainagePhaseRelPermMinEndPoint,
         drainagePhaseRelPermMaxEndPoint;
  // Step 1: Read the drainage for the non wetting phase
  GEOSX_ASSERT( m_drainageWettingNonWettingRelPermTableNames.size() == 2 );
  auto tableName = ( numPhases == 2 ) ? m_drainageWettingNonWettingRelPermTableNames[1] :
                   m_drainageNonWettingIntermediateRelPermTableNames[0];
  checkExistenceAndValidateRelPermTable( tableName,       // input
                                         drainagePhaseMinVolFraction,      // output
                                         drainagePhaseMaxVolFraction,
                                         drainagePhaseRelPermMinEndPoint,
                                         drainagePhaseRelPermMaxEndPoint );

  // Step 2: validate non-wetting-phase imbibition relative permeability table
  real64 imbibitionPhaseMinVolFraction, // output
         imbibitionPhaseMaxVolFraction,
         imbibitionPhaseRelPermMinEndPoint,
         imbibitionPhaseRelPermMaxEndPoint;

  if( m_phaseHasHysteresis[IPT::NONWETTING] )
  {

    checkExistenceAndValidateRelPermTable( m_imbibitionNonWettingRelPermTableName,  // input
                                           imbibitionPhaseMinVolFraction,  // output
                                           imbibitionPhaseMaxVolFraction,
                                           imbibitionPhaseRelPermMinEndPoint,
                                           imbibitionPhaseRelPermMaxEndPoint );

    GEOSX_THROW_IF( !isZero ( imbibitionPhaseMaxVolFraction - drainagePhaseMaxVolFraction ),
                    GEOSX_FMT( string( "{}: the maximum non-wetting-phase volume fraction (saturation) must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage maximum wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition maximum wetting-phase volume fraction is {}" ),
                               getFullName(),
                               drainagePhaseMaxVolFraction, imbibitionPhaseMaxVolFraction ),
                    InputError );

    GEOSX_THROW_IF( !isZero ( imbibitionPhaseRelPermMaxEndPoint - drainagePhaseRelPermMaxEndPoint ),
                    GEOSX_FMT( string( "{}: the non-wetting-phase relperm endpoint must be the same in drainage and imbibition.\n" )
                               + string( "However, we found that the drainage endpoint wetting-phase relperm is {}, " )
                               + string( "whereas the imbibition endpoint wetting-phase relperm is {}" ),
                               getFullName(),
                               drainagePhaseRelPermMaxEndPoint, imbibitionPhaseRelPermMaxEndPoint ),
                    InputError );

    GEOSX_THROW_IF( imbibitionPhaseMinVolFraction < drainagePhaseMinVolFraction,
                    GEOSX_FMT( string( "{}: the critical wetting-phase volume fraction (saturation) must be larger in imbibition (compared to the drainage value).\n" )
                               + string( "However, we found that the drainage critical wetting-phase volume fraction is {}, " )
                               + string( "whereas the imbibition critical wetting-phase volume fraction is {}" ),
                               getFullName(),
                               drainagePhaseMinVolFraction, imbibitionPhaseMinVolFraction ),
                    InputError );

    // test the opposite End point ?

  }


  m_nonWettingCurve.setPoints( {drainagePhaseMaxVolFraction, drainagePhaseRelPermMaxEndPoint},   // same as imbibition max
                               {imbibitionPhaseMinVolFraction, imbibitionPhaseRelPermMinEndPoint},
                               {drainagePhaseMinVolFraction, drainagePhaseRelPermMinEndPoint} );

}

void TableRelativePermeabilityHysteresis::checkExistenceAndValidateIntermediateRelPermTables()
{

  if( m_phaseNames.size() == 3 )
  {

    real64 drainagePhaseMinVolFraction,     //
           drainagePhaseMaxVolFraction,
           drainagePhaseRelPermMinEndPoint,
           drainagePhaseRelPermMaxEndPoint;


    //intermediate drainage from wetting
    checkExistenceAndValidateRelPermTable( m_drainageWettingIntermediateRelPermTableNames[1],    // input
                                           drainagePhaseMinVolFraction,    // output
                                           drainagePhaseMaxVolFraction,
                                           drainagePhaseRelPermMinEndPoint,
                                           drainagePhaseRelPermMaxEndPoint );

    //??
    checkExistenceAndValidateRelPermTable( m_drainageNonWettingIntermediateRelPermTableNames[1],    // input
                                           drainagePhaseMinVolFraction,
                                           drainagePhaseMaxVolFraction,
                                           drainagePhaseRelPermMinEndPoint,
                                           drainagePhaseRelPermMaxEndPoint );
  }
}


void TableRelativePermeabilityHysteresis::checkExistenceAndValidateRelPermTable( string const & relPermTableName,
                                                                                 real64 & phaseMinVolFrac,
                                                                                 real64 & phaseMaxVolFrac,
                                                                                 real64 & phaseRelPermMinEndPoint,
                                                                                 real64 & phaseRelPermMaxEndPoint ) const
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
  TableRelativePermeabilityHelpers::validateRelativePermeabilityTable( relPermTable,  // input
                                                                       fullName,
                                                                       phaseMinVolFrac,  // output
                                                                       phaseMaxVolFrac,
                                                                       phaseRelPermMinEndPoint,
                                                                       phaseRelPermMaxEndPoint );
}

void TableRelativePermeabilityHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;
  m_KilloughModel.computeLandCoefficient( m_wettingCurve, m_landParam[IPT::WETTING] );
  m_KilloughModel.computeLandCoefficient( m_nonWettingCurve, m_landParam[IPT::NONWETTING] );
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

void TableRelativePermeabilityHysteresis::createKilloughKernelWrapper()
{

  m_KilloughKernel = m_KilloughModel.createKernelWrapper( m_landParam,
                                                          m_phaseTrappedVolFrac );

}

TableRelativePermeabilityHysteresis::KernelWrapper
TableRelativePermeabilityHysteresis::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  createKilloughKernelWrapper();

  // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
  return KernelWrapper( m_drainageRelPermKernelWrappers,
                        m_imbibitionRelPermKernelWrappers,
                        m_KilloughKernel,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_wettingCurve,
                        m_nonWettingCurve,
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
                                                                   KilloughHysteresis::KernelKilloughHysteresisBase const & killoughHystWrappers,
                                                                   arrayView1d< integer const > const & phaseHasHysteresis,
                                                                   arrayView1d< real64 const > const & landParam,
                                                                   KilloughHysteresis::HysteresisCurve const & wettingCurve,
                                                                   KilloughHysteresis::HysteresisCurve const & nonWettingCurve,
                                                                   arrayView1d< integer const > const & phaseTypes,
                                                                   arrayView1d< integer const > const & phaseOrder,
                                                                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
                                                                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
                                                                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrappedVolFrac,
                                                                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                                                   arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
  :
  RelativePermeabilityBaseUpdate( phaseTypes,
                                  phaseOrder,
                                  phaseRelPerm,
                                  dPhaseRelPerm_dPhaseVolFrac,
                                  phaseTrappedVolFrac ),
  m_KilloughKernel( killoughHystWrappers ),
  m_drainageRelPermKernelWrappers( drainageRelPermKernelWrappers ),
  m_imbibitionRelPermKernelWrappers( imbibitionRelPermKernelWrappers ),
  m_phaseHasHysteresis( phaseHasHysteresis ),
  m_landParam( landParam ),
  m_wettingCurve( wettingCurve ),
  m_nonWettingCurve( nonWettingCurve ),
  m_phaseMinHistoricalVolFraction( phaseMinHistoricalVolFraction ),
  m_phaseMaxHistoricalVolFraction( phaseMaxHistoricalVolFraction )
{}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeabilityHysteresis, std::string const &, Group * const )

} // namespace constitutive

} // namespace geosx
