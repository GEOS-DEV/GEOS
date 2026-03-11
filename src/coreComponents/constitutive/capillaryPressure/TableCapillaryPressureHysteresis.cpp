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


#include "TableCapillaryPressureHysteresis.hpp"

#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

TableCapillaryPressureHysteresis::TableCapillaryPressureHysteresis( const std::string & name,
                                                                    dataRepository::Group * const parent )
  : CapillaryPressureBase( name, parent )
{

  registerWrapper( viewKeyStruct::phaseHasHysteresisString(), &m_phaseHasHysteresis ).
    setInputFlag( InputFlags::FALSE )
    .                                         // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::landParameterString(), &m_landParam ).
    setInputFlag( InputFlags::FALSE ).                                     // will be deduced from tables
    setSizedFromParent( 0 );

  //2phase
  registerWrapper( viewKeyStruct::drainageWettingNonWettingCapPresTableNameString(),
                   &m_drainageWettingNonWettingCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the drainage two-phase table for capillary pressure curve. \n"
                    "If you want to use 3-phase flow please use instead " +
                    string( viewKeyStruct::drainageWettingIntermediateCapPresTableNameString()) +
                    " and " +
                    string( viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString()) +
                    "to specify the tables names" );
  registerWrapper( viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(),
                   &m_imbibitionWettingNonWettingCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the drainage two-phase table for capillary pressure curve. \n"
                    "If you want to use 3-phase flow please use instead " +
                    string( viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString()) +
                    " and " +
                    string( viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString()) +
                    "to specify the tables names" );
  //3phase
  registerWrapper( viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),
                   &m_drainageWettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Drainage wetting/intermediate (e.g. w/o) capillary pressure table name for the wetting phase.\n"
    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );
  registerWrapper( viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString(),
                   &m_drainageNonWettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Drainage non-wetting/intermediate (e.g. o/g) capillary pressure table name for the non-wetting phase.\n"
    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );
  registerWrapper( viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(),
                   &m_imbibitionWettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Imbibition wetting/intermediate (e.g. w/o) table name for the wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );
  registerWrapper( viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString(),
                   &m_imbibitionNonWettingIntermediateCapPresTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Imbibition non-wetting/intermediate (e.g. o/g) table name for the wetting phase.\n"
                    "To neglect hysteresis on this phase, just use the same table name for the drainage and imbibition curves" );

  // kernels
  //2p
  registerWrapper( viewKeyStruct::wettingNonWettingCapillaryPressureKernelWrappersString(),
                   &m_wettingNonWettingCapillaryPressureKernelWrappers )
    .setSizedFromParent( 0 ).setRestartFlags( RestartFlags::NO_WRITE );
  //3p
  registerWrapper( viewKeyStruct::wettingIntermediateCapillaryPressureKernelWrappersString(),
                   &m_wettingIntermediateCapillaryPressureKernelWrappers )
    .setSizedFromParent( 0 ).setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::nonWettingIntermediateCapillaryPressureKernelWrappersString(),
                   &m_nonWettingIntermediateCapillaryPressureKernelWrappers )
    .setSizedFromParent( 0 ).setRestartFlags( RestartFlags::NO_WRITE );


  registerWrapper( viewKeyStruct::wettingCurveString(), &m_wettingCurve ).
    setInputFlag(
    InputFlags::FALSE ).                                // will be deduced from tables
    setSizedFromParent(
    0 )
    .setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::nonWettingCurveString(), &m_nonWettingCurve ).
    setInputFlag(
    InputFlags::FALSE ).                                // will be deduced from tables
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


  registerWrapper( KilloughHysteresis::viewKeyStruct::killoughCurvatureParameterPcString(),
                   &m_killoughCurvatureParamCapPres ).
    setInputFlag(
    InputFlags::OPTIONAL ).
    setApplyDefaultValue(
    .1 ).
    setDescription(
    "Curvature parameter introduced by Killough for wetting-phase hysteresis (see RTD documentation)." );

  //misc
  registerWrapper( viewKeyStruct::phaseIntermediateMinVolFractionString(), &m_phaseIntermediateMinVolFraction ).
    setInputFlag( InputFlags::FALSE ).setDescription( "min vol fraction of intermediate if exist" ).
    // will be deduced from tables
    setSizedFromParent( 0 );

  registerField< fields::cappres::mode >( &m_mode );


  registerField< fields::cappres::phaseMaxHistoricalVolFraction >(
    &m_phaseMaxHistoricalVolFraction );
  registerField< fields::cappres::phaseMinHistoricalVolFraction >(
    &m_phaseMinHistoricalVolFraction );
  registerField< fields::cappres::phaseMode2PeakVolFraction >(
    &m_phaseMode2PeakVolFraction );

}

/// usual utils

void TableCapillaryPressureHysteresis::postProcessInput()
{

  using TPP = ThreePhasePairPhaseType;

  integer const numPhases = m_phaseNames.size();
  GEOS_THROW_IF( numPhases != 2 && numPhases != 3,
                 GEOS_FMT( "{}: the expected number of fluid phases is either two, or three",
                           getFullName()),
                 InputError );

  m_phaseHasHysteresis.resize( 2 );

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( m_drainageWettingNonWettingCapPresTableName.empty(),
                   GEOS_FMT(
                     "{}: for a two-phase flow simulation, we must use {} to specify the capillary pressure table for the drainage pair (wetting phase, non-wetting phase)",
                     getFullName(),
                     viewKeyStruct::drainageWettingNonWettingCapPresTableNameString()),
                   InputError );


    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = (m_imbibitionWettingNonWettingCapPresTableName.empty() ||
                                                       m_imbibitionWettingNonWettingCapPresTableName ==
                                                       m_drainageWettingNonWettingCapPresTableName)
                                                                  ? 0 : 1;
    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING];


  }
  else if( numPhases == 3 )
  {


    GEOS_THROW_IF( m_drainageWettingIntermediateCapPresTableName.empty() ||
                   m_drainageNonWettingIntermediateCapPresTableName.empty(),
                   GEOS_FMT(
                     "{}: for a three-phase flow simulation, we must use {} to specify the capillary pressure table "
                     "for the pair (wetting phase, intermediate phase), and {} to specify the capillary pressure table "
                     "for the pair (non-wetting phase, intermediate phase)",
                     getFullName(),
                     viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),
                     viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString()),
                   InputError );

    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = (m_imbibitionWettingIntermediateCapPresTableName.empty() ||
                                                       m_imbibitionWettingIntermediateCapPresTableName ==
                                                       m_drainageWettingIntermediateCapPresTableName)
                                                                  ? 0 : 1;

    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = (m_imbibitionNonWettingIntermediateCapPresTableName.empty() ||
                                                          m_imbibitionNonWettingIntermediateCapPresTableName ==
                                                          m_drainageNonWettingIntermediateCapPresTableName)
                                                                     ? 0 : 1;
  }
  //Killough section
  //TODO improve hard coded default
  KilloughHysteresis::postProcessInput( m_jerauldParam_a, m_jerauldParam_b, 0,
                                        m_killoughCurvatureParamCapPres );

  GEOS_THROW_IF( m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] == 0 &&
                 m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] == 0,
                 GEOS_FMT(
                   "{}: we must use {} (2-phase) / {} or {} (3-phase) to specify at least one imbibition relative permeability table",
                   getFullName(),
                   viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(),
                   viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(),
                   viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString()),
                 InputError );

}

void TableCapillaryPressureHysteresis::initializePreSubGroups()
{
  CapillaryPressureBase::initializePreSubGroups();

  integer const numPhases = m_phaseNames.size();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  //equivalent to oil/gas - a.k.a two phase flow ordered by non wetting
  bool const capPresMustBeIncreasing = (m_phaseOrder[PhaseType::WATER] < 0)
                                                 ? true   // pc on the gas phase, function must be increasing
                                                 : false; // pc on the water phase, function must be decreasing


  // Step 1: check sanity of drainage tables
  if( numPhases == 2 )
  {

    real64 drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres,
           drainageNonWettingPhaseMinVolumeFraction, drainageNonWettingMinCapPres,
           imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres,
           imbibitionNonWettingPhaseMinVolumeFraction, imbibitionNonWettingMinCapPres,
           wettingPhaseMinVolumeFraction, wettingMaxCapPres,
           nonWettingPhaseMaxVolumeFraction, nonWettingMaxCapPres;

    {

      imbibitionNonWettingMinCapPres = 0.0;

      GEOS_THROW_IF( !functionManager.hasGroup( m_drainageWettingNonWettingCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_drainageWettingNonWettingCapPresTableName ),
                     InputError );
      TableFunction const
      & capPresTable = functionManager.getGroup< TableFunction >(
        m_drainageWettingNonWettingCapPresTableName );

      //w/o  or  w/g pair
      if( !capPresMustBeIncreasing )
      {
        TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                       capPresMustBeIncreasing,
                                                                       drainageWettingPhaseMaxVolumeFraction,
                                                                       wettingPhaseMinVolumeFraction,
                                                                       drainageWettingMinCapPres,
                                                                       wettingMaxCapPres );

        drainageNonWettingPhaseMinVolumeFraction = 1. - drainageWettingPhaseMaxVolumeFraction;
        nonWettingPhaseMaxVolumeFraction = 1. - wettingPhaseMinVolumeFraction;

      }
      else                   // o/g pair
      {
        TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                       capPresMustBeIncreasing,
                                                                       nonWettingPhaseMaxVolumeFraction,
                                                                       drainageNonWettingPhaseMinVolumeFraction,
                                                                       nonWettingMaxCapPres,
                                                                       drainageNonWettingMinCapPres );

        drainageWettingPhaseMaxVolumeFraction = 1. - drainageNonWettingPhaseMinVolumeFraction;
        wettingPhaseMinVolumeFraction = 1. - nonWettingPhaseMaxVolumeFraction;
      }

    }

    {
      GEOS_THROW_IF( !functionManager.hasGroup( m_imbibitionWettingNonWettingCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_imbibitionWettingNonWettingCapPresTableName ),
                     InputError );
      TableFunction const
      & capPresTable = functionManager.getGroup< TableFunction >(
        m_imbibitionWettingNonWettingCapPresTableName );

      //w/o  or  w/g pair
      if( !capPresMustBeIncreasing )
      {
        TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                       capPresMustBeIncreasing,
                                                                       imbibitionWettingPhaseMaxVolumeFraction,
                                                                       wettingPhaseMinVolumeFraction,
                                                                       imbibitionWettingMinCapPres,
                                                                       wettingMaxCapPres );

        imbibitionNonWettingPhaseMinVolumeFraction = 1. - imbibitionWettingPhaseMaxVolumeFraction;
        nonWettingPhaseMaxVolumeFraction = 1. - wettingPhaseMinVolumeFraction;

      }
      else                   // o/g pair
      {
        TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                       capPresMustBeIncreasing,
                                                                       nonWettingPhaseMaxVolumeFraction,
                                                                       imbibitionNonWettingPhaseMinVolumeFraction,
                                                                       nonWettingMaxCapPres,
                                                                       imbibitionWettingMinCapPres );

        imbibitionWettingPhaseMaxVolumeFraction = 1. - imbibitionNonWettingPhaseMinVolumeFraction;
        wettingPhaseMinVolumeFraction = 1. - nonWettingPhaseMaxVolumeFraction;
      }
    }

    //constructing wetting/nonwetting curves

    if( !capPresMustBeIncreasing )
    {
      m_wettingCurve.setPoints(
        {wettingPhaseMinVolumeFraction, wettingMaxCapPres},                       // same as imbibition min
        {imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres},
        {drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres} );
    }
    else
    {
      m_nonWettingCurve.setPoints(
        {nonWettingPhaseMaxVolumeFraction, nonWettingMaxCapPres},
        {imbibitionNonWettingPhaseMinVolumeFraction, imbibitionNonWettingMinCapPres},
        {drainageNonWettingPhaseMinVolumeFraction, drainageNonWettingMinCapPres}
        );
    }

  }
  else if( numPhases == 3 )
  {

    real64 drainageWettingPhaseMaxVolumeFraction, drainageWettingMinCapPres,
           drainageNonWettingPhaseMinVolumeFraction, drainageNonWettingMinCapPres,
           imbibitionWettingPhaseMaxVolumeFraction, imbibitionWettingMinCapPres,
           imbibitionNonWettingPhaseMinVolumeFraction, imbibitionNonWettingMinCapPres,
           wettingPhaseMinVolumeFraction, wettingMaxCapPres,
           nonWettingPhaseMaxVolumeFraction, nonWettingMaxCapPres;

    GEOS_UNUSED_VAR( drainageWettingMinCapPres );
//define scope to avoid differentiate temp var (lazy)
    {
      GEOS_THROW_IF( !functionManager.hasGroup( m_drainageWettingIntermediateCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_drainageWettingIntermediateCapPresTableName ),
                     InputError );
      TableFunction const
      & capPresTableWI = functionManager.getGroup< TableFunction >(
        m_drainageWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableWI, getFullName(), false,
                                                                     drainageWettingPhaseMaxVolumeFraction,
                                                                     wettingPhaseMinVolumeFraction,
                                                                     drainageNonWettingMinCapPres,
                                                                     wettingMaxCapPres );

      GEOS_THROW_IF( !functionManager.hasGroup( m_drainageNonWettingIntermediateCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_drainageNonWettingIntermediateCapPresTableName ),
                     InputError );
      TableFunction const & capPresTableNWI =
        functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableNWI, getFullName(), true,
                                                                     nonWettingPhaseMaxVolumeFraction,
                                                                     drainageNonWettingPhaseMinVolumeFraction,
                                                                     nonWettingMaxCapPres,
                                                                     drainageWettingPhaseMaxVolumeFraction
                                                                     );

      m_phaseIntermediateMinVolFraction =
        1.0 - drainageWettingPhaseMaxVolumeFraction - drainageWettingPhaseMaxVolumeFraction;
    }

    if( !m_imbibitionWettingIntermediateCapPresTableName.empty())
    {

      GEOS_THROW_IF( !functionManager.hasGroup( m_imbibitionWettingIntermediateCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_imbibitionWettingIntermediateCapPresTableName ),
                     InputError );
      TableFunction const
      & capPresTableWI = functionManager.getGroup< TableFunction >(
        m_imbibitionWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableWI, getFullName(), false,
                                                                     imbibitionWettingPhaseMaxVolumeFraction,
                                                                     wettingPhaseMinVolumeFraction,
                                                                     imbibitionWettingMinCapPres,
                                                                     wettingMaxCapPres
                                                                     );


    }

    if( !m_imbibitionNonWettingIntermediateCapPresTableName.empty())
    {

      GEOS_THROW_IF( !functionManager.hasGroup( m_imbibitionNonWettingIntermediateCapPresTableName ),
                     GEOS_FMT( "{}: the table function named {} could not be found",
                               getFullName(),
                               m_imbibitionNonWettingIntermediateCapPresTableName ),
                     InputError );
      TableFunction const & capPresTableNWI =
        functionManager.getGroup< TableFunction >( m_imbibitionNonWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableNWI, getFullName(), true,
                                                                     nonWettingPhaseMaxVolumeFraction,
                                                                     imbibitionNonWettingPhaseMinVolumeFraction,
                                                                     nonWettingMaxCapPres,
                                                                     imbibitionNonWettingMinCapPres );


    }
  }

  // Step 2: check the sanity btw drainage and imbibition
  auto const eps = 1e-15;
  if( numPhases == 2 )
  {
    //TODO weak make stronger
    GEOS_THROW_IF(
      m_wettingCurve.isZero() && m_nonWettingCurve.isZero(),
      GEOS_FMT(
        "{}: Inconsistent data for capillary pressure hysteresis. No hysteresis curve is defined.",
        getFullName()),
      InputError );

    GEOS_THROW_IF(
      !m_wettingCurve.isZero() && !m_nonWettingCurve.isZero(),
      GEOS_FMT(
        "{}: Inconsistent data for capillary pressure hysteresis. Both non wetting and wetting hysteresis curve are defined in two phase flow setting.",
        getFullName()),
      InputError );


  }
  else if( numPhases == 3 )
  {

    GEOS_THROW_IF( std::fabs( m_wettingCurve.oppositeBoundPhaseVolFraction - (1. - m_nonWettingCurve.oppositeBoundPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                   GEOS_FMT(
                     "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                     getFullName(), "Sw_min", "Snw_max", "Sinter_min" ),
                   InputError );
    GEOS_THROW_IF( std::fabs( m_wettingCurve.drainageExtremaPhaseVolFraction - (1. - m_nonWettingCurve.drainageExtremaPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                   GEOS_FMT(
                     "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                     getFullName(), "Sw_min", "Snw_max", "Sinter_min" ),
                   InputError );
    GEOS_THROW_IF( std::fabs( m_wettingCurve.imbibitionExtremaPhaseVolFraction - (1. - m_nonWettingCurve.imbibitionExtremaPhaseVolFraction - m_phaseIntermediateMinVolFraction)) > eps,
                   GEOS_FMT(
                     "{}: Inconsistent data for capillary pressure hysteresis. {}, {} and {} should sum up to 1.",
                     getFullName(), "Sw_min", "Snw_max", "Sinter_min" ),
                   InputError );

  }


  // Step 3: compute the Land coefficient
  computeLandCoefficient();



  if( m_phaseMaxHistoricalVolFraction.size( 1 ) == 0 && numPhases > 0 )
  {
    localIndex const currentSize = m_phaseMaxHistoricalVolFraction.size( 0 );
    if( currentSize > 0 )
    {
      m_phaseMaxHistoricalVolFraction.resize( currentSize, numPhases );
      m_phaseMinHistoricalVolFraction.resize( currentSize, numPhases );
      m_phaseMode2PeakVolFraction.resize( currentSize, numPhases );
      m_phaseMaxHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
      m_phaseMinHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 1.0 );
      m_phaseMode2PeakVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
    }
  }
}

/// Land coeff (tb refactored out in KilloughHysteresis) and saved cvgd

void TableCapillaryPressureHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase

  integer ipWetting, ipNonWetting;
  std::tie( ipWetting, ipNonWetting ) = phaseIndex( m_phaseOrder );

  KilloughHysteresis::computeLandCoefficient( m_wettingCurve, m_landParam[ipWetting] );
  KilloughHysteresis::computeLandCoefficient( m_nonWettingCurve, m_landParam[ipNonWetting] );

}

/// common utils
void TableCapillaryPressureHysteresis::resizeFields( localIndex const size, localIndex const numPts )
{
  CapillaryPressureBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();



  m_mode.resize( size );

  if( numPhases > 0 )
  {
    m_phaseMaxHistoricalVolFraction.resize( size, numPhases );
    m_phaseMinHistoricalVolFraction.resize( size, numPhases );
    m_phaseMode2PeakVolFraction.resize( size, numPhases );
    m_phaseMaxHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
    m_phaseMinHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 1.0 );
    m_phaseMode2PeakVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
  }


}

void TableCapillaryPressureHysteresis::saveConvergedPhaseVolFractionState(
  arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
{
  CapillaryPressureBase::saveConvergedState();

  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMode2PeakVolFraction = m_phaseMode2PeakVolFraction.toView();
  arrayView1d< integer > mode_int = m_mode.toView();

  localIndex const numElems = phaseVolFraction.size( 0 );
  integer const numPhases = numFluidPhases();



  using PT = CapillaryPressureBase::PhaseType;
  arrayView1d< integer const > const phaseOrderView = phaseOrder();
  integer const ipWater = phaseOrderView[PT::WATER];
  integer const ipOil = phaseOrderView[PT::OIL];
  integer const ipGas = phaseOrderView[PT::GAS];

  integer ipWetting = -1;
  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    ipWetting = ipWater;
  }
  else if( ipWater < 0 )
  {
    ipWetting = ipOil;
  }
  else if( ipOil < 0 )
  {
    ipWetting = ipWater;
  }
  else if( ipGas < 0 )
  {
    ipWetting = ipWater;
  }

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei ) {
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMaxHistoricalVolFraction[ei][ip] = LvArray::math::max( phaseVolFraction[ei][ip],
                                                                  phaseMaxHistoricalVolFraction[ei][ip] );
      phaseMinHistoricalVolFraction[ei][ip] = LvArray::math::min( phaseVolFraction[ei][ip],
                                                                  phaseMinHistoricalVolFraction[ei][ip] );
    }

    if( ipWetting >= 0 && ipWetting < numPhases )
    {
      ModeIndexType const currentMode = static_cast< ModeIndexType >(mode_int[ei]);
      if( currentMode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
      {
        real64 const currentS = phaseVolFraction[ei][ipWetting];
        if( currentS > phaseMode2PeakVolFraction[ei][ipWetting] )
        {
          phaseMode2PeakVolFraction[ei][ipWetting] = currentS;
        }
      }
      else if( currentMode != ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
      {
        phaseMode2PeakVolFraction[ei][ipWetting] = 0.0;
      }
    }
  } );

}

void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionWettingCapillaryPressure(
  const arrayView1d< const TableFunction::KernelWrapper > & wettingKernelWapper,
  const KilloughHysteresis::HysteresisCurve & wettingCurve,
  const KilloughHysteresis::HysteresisCurve & nonWettingCurve,              //discard if not needed
  const geos::real64 & landParam,
  const geos::real64 & phaseVolFraction,
  const geos::real64 & phaseMinHistoricalVolFraction,
  const geos::real64 & phaseMaxHistoricalVolFraction,
  const geos::real64 & phaseMode2PeakVolFraction,
  geos::real64 & phaseTrappedVolFrac,
  geos::real64 & phaseCapPressure,
  geos::real64 & dPhaseCapPressure_dPhaseVolFrac,
  const ModeIndexType & mode ) const
{
  GEOS_ASSERT( wettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smxi = wettingCurve.imbibitionExtremaPhaseVolFraction;
  real64 const Smxd = wettingCurve.drainageExtremaPhaseVolFraction;
  real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;
  real64 const Smax = wettingCurve.drainageExtremaPhaseVolFraction;

  GEOS_UNUSED_VAR( Smxi, Smxd, phaseTrappedVolFrac );

//  if( S <= Smin )
//  {
//    //below accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else if( S >= Smxd )
//  {
//    //above accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = wettingKernelWapper[ModeIndexType::IMBIBITION].compute( &S, &dpci_dS );
    real64 const pcd = wettingKernelWapper[ModeIndexType::DRAINAGE].compute( &S, &dpcd_dS );
    real64 const Somin = m_phaseIntermediateMinVolFraction;

    real64 const E = m_killoughCurvatureParamCapPres;

    //Step 2. compute F as in (EQ 34.15) F = (1/(Sw-Shy+E)-1/E) / (1/(Swma-Shy+E)-1/E)
    //drainage to imbibition branch
    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
    {
      // DRAINAGE_TO_IMBIBITION: Shy should be minimum historical (where drainage started)
      real64 const Shy = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;

      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( nonWettingCurve,
                                                                  Shy,
                                                                  landParam,
                                                                  m_jerauldParam_a,
                                                                  m_jerauldParam_b,
                                                                  Scrt );
      real64 const Swma = 1 - Scrt - Somin;
      real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swma - Shy + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. Eventually assemble everything following (EQ. 34.14)
      phaseCapPressure = pcd + F * (pci - pcd);
      dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpci_dS - dpcd_dS);
    }
    //imbibition to drainage
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      // IMBIBITION_TO_DRAINAGE: Shy should be maximum historical (where imbibition started)
      real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

      // For IMBIBITION_TO_DRAINAGE, use the same formula structure as non-wetting phase
      // F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swmin - Shy + E) - 1. / E)
      // This ensures F = 0 when S = Shy (high saturation) and F = 1 when S = Swmin (low saturation)
      // Minimum accessible wetting saturation = Somin (irreducible wetting saturation)
      // Use actual Somin (may be 0.0), the formula will handle negative values correctly
      real64 const Swmin = Somin;

      real64 const F_num = (1. / (S - Shy + E) - 1. / E);
      real64 const F_denom = (1. / (Swmin - Shy + E) - 1. / E);
      // Both numerator and denominator can be negative when Swmin < Shy - E
      // This is okay - negative/negative = positive F
      real64 F = F_num / F_denom;
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. Eventually assemble everything following (EQ. 34.14)
      phaseCapPressure = pci + F * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpcd_dS - dpci_dS);
    }
    //imbibition to drainage from scanning curve (departing from DRAINAGE_TO_IMBIBITION)
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
    {
      // IMBIBITION_TO_DRAINAGE_FROM_SCANNING: Secondary drainage scanning curve
      // Based on Killough (1976) - solve quadratic for ghost departure point
      //
      // Nomenclature:
      // - Sw_star (S_star): second reversal saturation (turnaround point, peak reached during Mode 2)
      // - Sw_Hyst_imb (H): first reversal point (where imbibition started, departure from drainage)
      // - Sw_wr (Sw_wr): residual/connate wetting saturation
      // - x (Sw_star_Hyst_ghost): ghost departure point on imbibition curve (unknown, to be solved)
      //
      // Step 1: Identify reversal points
      // S_star is the turnaround point (peak reached during Mode 2)
      real64 S_star = phaseMode2PeakVolFraction;
      if( S_star < 1e-12 || S_star >= Smax - 1e-12 )
      {
        // Use phaseMaxHistoricalVolFraction, but ensure it's between H and Smax
        real64 const H_temp = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;
        if( phaseMaxHistoricalVolFraction > H_temp + 1e-12 && phaseMaxHistoricalVolFraction < Smax - 1e-12 )
        {
          S_star = phaseMaxHistoricalVolFraction;
        }
        else
        {
          S_star = LvArray::math::max( H_temp + 0.1, LvArray::math::min( S, Smax - 0.01 ));
        }
      }

      real64 const H = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;
      // Sw_wr: residual/connate wetting saturation - minimum saturation on drainage curve
      // For wetting phase, oppositeBoundPhaseVolFraction should be the minimum, but if it's 0,
      // we need to get the actual minimum from the drainage curve table
      // Since we can't easily access table coordinates from KernelWrapper, use Smin if > 0, otherwise try Somin
      real64 Sw_wr = Smin;
      if( Sw_wr < 1e-12 )
      {
        Sw_wr = (Somin > 1e-12) ? Somin : 1e-6;                 // Small default to avoid division issues
      }

      // Step 2: Compute F_known from Mode 2 (imbibition scanning curve) at S_star
      // Mode 2 formula: F = (1 / (S - H + E) - 1 / E) / (1 / (Swma - H + E) - 1 / E)
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( nonWettingCurve,
                                                                  H,
                                                                  landParam,
                                                                  m_jerauldParam_a,
                                                                  m_jerauldParam_b,
                                                                  Scrt );
      real64 const Swma = 1 - Scrt - Somin;
      real64 F_known = 0.0;
      if( S_star <= Swma + 1e-12 )
      {
        real64 const F_num = (1. / (S_star - H + E) - 1. / E);
        real64 const F_denom = (1. / (Swma - H + E) - 1. / E);
        if( LvArray::math::abs( F_denom ) > 1e-12 )
        {
          F_known = F_num / F_denom;
        }
        else
        {
          F_known = (S_star > H) ? 1.0 : 0.0;
        }
      }
      else
      {
        // S_star > Swma: beyond Mode 2 range, use F_known = 0 (pure imbibition)
        F_known = 0.0;
      }
      F_known = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_known ));

      // Step 3: Compute F_star_target = 1 - F_known (from Eq. A-6)
      real64 const F_star_target = 1.0 - F_known;

      // Step 4: Solve quadratic for x = Sw_star_Hyst_ghost
      // F* = [1/(x - S_star + E) - 1/E] / [1/(x - H + E) - 1/E] = F_star_target
      // Note: Denominator uses H (first reversal point), not Sw_wr, so curve rejoins drainage at H
      // Rearranging gives: (1-T)*x^2 + B*x + C = 0
      // where T = F_star_target, a = S_star, b = H, e = E
      real64 const a = S_star;
      real64 const b = H;                // Use H instead of Sw_wr so curve rejoins drainage at H
      real64 const e = E;
      real64 const T = F_star_target;

      real64 const coeff_A = (1.0 - T);
      real64 const coeff_B = -(1.0 - T) * (a + b - 2.0 * e) - (1.0 - T) * e;
      real64 const coeff_C = (1.0 - T) * (a - e) * (b - e)
                             - e * e * (1.0 + T)
                             - e * (T * a - b);

      real64 x = S_star;
      real64 discriminant = coeff_B * coeff_B - 4.0 * coeff_A * coeff_C;

      if( discriminant >= 0.0 && LvArray::math::abs( coeff_A ) > 1e-14 )
      {
        real64 sqrt_disc = LvArray::math::sqrt( discriminant );
        real64 x1 = (-coeff_B + sqrt_disc) / (2.0 * coeff_A);
        real64 x2 = (-coeff_B - sqrt_disc) / (2.0 * coeff_A);

        // Pick the physically meaningful root: x must be > S_star
        if( x1 > S_star - 1e-8 && x2 > S_star - 1e-8 )
        {
          x = LvArray::math::min( x1, x2 );                 // Take the closer one
        }
        else if( x1 > S_star - 1e-8 )
        {
          x = x1;
        }
        else if( x2 > S_star - 1e-8 )
        {
          x = x2;
        }
      }

      // Step 5: Compute F_star for current saturation S
      // F_star = [1/(x - S + E) - 1/E] / [1/(x - H + E) - 1/E]
      // When S = H, F_star = 1.0 (pure drainage), so curve rejoins drainage at H
      real64 F_star = 1.0;
      real64 F_star_num = (1. / (x - S + E) - 1. / E);
      real64 F_star_denom = (1. / (x - H + E) - 1. / E);
      if( LvArray::math::abs( F_star_denom ) > 1e-12 )
      {
        F_star = F_star_num / F_star_denom;
        F_star = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_star ));
      }

      // Step 6: Compute Pc using Eq. A-5: P_c = P_c^Im + F* [P_c^Dr - P_c^Im]
      phaseCapPressure = pci + F_star * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F_star * (dpcd_dS - dpci_dS);
    }
    else
    {
      GEOS_THROW( GEOS_FMT( "{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                            "TableCapillaryPressureHysteresis",
                            (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                               ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                  InputError );
    }


  }

}

void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionWettingCapillaryPressure(
  const arrayView1d< const TableFunction::KernelWrapper > & wettingKernelWapper,
  const KilloughHysteresis::HysteresisCurve & wettingCurve,
  const geos::real64 & landParam,
  const geos::real64 & phaseVolFraction,
  const geos::real64 & phaseMinHistoricalVolFraction,
  const geos::real64 & phaseMaxHistoricalVolFraction,
  const geos::real64 & phaseMode2PeakVolFraction,
  geos::real64 & phaseTrappedVolFrac,
  geos::real64 & phaseCapPressure,
  geos::real64 & dPhaseCapPressure_dPhaseVolFrac,
  const ModeIndexType & mode ) const
{
  GEOS_ASSERT( wettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smxi = wettingCurve.imbibitionExtremaPhaseVolFraction;
  real64 const Smxd = wettingCurve.drainageExtremaPhaseVolFraction;
  real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;
  real64 const Smax = wettingCurve.drainageExtremaPhaseVolFraction;

  GEOS_UNUSED_VAR( Smxi, Smxd, phaseTrappedVolFrac );

//  if( S <= Smin )
//  {
//    //below accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else if( S >= Smxd )
//  {
//    //above accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = -CAP_INF_DERIV;
//  }
//  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = wettingKernelWapper[ModeIndexType::IMBIBITION].compute( &S, &dpci_dS );
    real64 const pcd = wettingKernelWapper[ModeIndexType::DRAINAGE].compute( &S, &dpcd_dS );

    real64 const Somin = m_phaseIntermediateMinVolFraction;
    real64 const E = m_killoughCurvatureParamCapPres;

    //Step 2. compute F as in (EQ 34.15) F = (1/(Sw-Shy+E)-1/E) / (1/(Swma-Shy+E)-1/E)
    //drainage to imbibition branch
    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
    {
      // DRAINAGE_TO_IMBIBITION: Shy should be minimum historical (where drainage started)
      real64 const Shy = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;

      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( wettingCurve,
                                                                  Shy,
                                                                  landParam,
                                                                  m_jerauldParam_a,
                                                                  m_jerauldParam_b,
                                                                  Scrt );



      //should be the pore space accessible to the two wetting phase
      real64 const Swma = 1 - (1 -  Scrt);
      real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swma - Shy + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. Eventually assemble everything following (EQ. 34.14)
      phaseCapPressure = pcd + F * (pci - pcd);
      dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpci_dS - dpcd_dS);
    }
    //imbibition to drainage
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      // IMBIBITION_TO_DRAINAGE: Shy should be maximum historical (where imbibition started)
      real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

      // For IMBIBITION_TO_DRAINAGE, use the same formula structure as non-wetting phase
      // F = (1. / (S - Shy + E) - 1. / E) / (1. / (Swmin - Shy + E) - 1. / E)
      // This ensures F = 0 when S = Shy (high saturation) and F = 1 when S = Swmin (low saturation)
      // For two-phase, use Smin as the minimum accessible wetting saturation
      // Use actual Smin (may be 0.0), the formula will handle negative values correctly
      real64 const Swmin = Smin;

      real64 const F_num = (1. / (S - Shy + E) - 1. / E);
      real64 const F_denom = (1. / (Swmin - Shy + E) - 1. / E);
      // Both numerator and denominator can be negative when Swmin < Shy - E
      // This is okay - negative/negative = positive F
      real64 F = F_num / F_denom;
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. Eventually assemble everything following (EQ. 34.14)
      phaseCapPressure = pci + F * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpcd_dS - dpci_dS);
    }
    //imbibition to drainage from scanning curve (departing from DRAINAGE_TO_IMBIBITION)
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
    {
      // IMBIBITION_TO_DRAINAGE_FROM_SCANNING: Secondary drainage scanning curve
      // Based on Killough (1976) - solve quadratic for ghost departure point
      //
      // Nomenclature:
      // - Sw_star (S_star): second reversal saturation (turnaround point, peak reached during Mode 2)
      // - Sw_Hyst_imb (H): first reversal point (where imbibition started, departure from drainage)
      // - Sw_wr (Sw_wr): residual/connate wetting saturation
      // - x (Sw_star_Hyst_ghost): ghost departure point on imbibition curve (unknown, to be solved)
      //
      // Step 1: Identify reversal points
      // S_star is the turnaround point (peak reached during Mode 2)
      real64 S_star = phaseMode2PeakVolFraction;
      if( S_star < 1e-12 || S_star >= Smax - 1e-12 )
      {
        // Use phaseMaxHistoricalVolFraction, but ensure it's between H and Smax
        real64 const H_temp = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;
        if( phaseMaxHistoricalVolFraction > H_temp + 1e-12 && phaseMaxHistoricalVolFraction < Smax - 1e-12 )
        {
          S_star = phaseMaxHistoricalVolFraction;
        }
        else
        {
          S_star = LvArray::math::max( H_temp + 0.1, LvArray::math::min( S, Smax - 0.01 ));
        }
      }

      real64 const H = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;
      // Sw_wr: residual/connate wetting saturation - minimum saturation on drainage curve
      // For wetting phase, oppositeBoundPhaseVolFraction should be the minimum, but if it's 0,
      // we need to get the actual minimum from the drainage curve table
      // Since we can't easily access table coordinates from KernelWrapper, use Smin if > 0, otherwise try Somin
      real64 Sw_wr = Smin;
      if( Sw_wr < 1e-12 )
      {
        Sw_wr = (Somin > 1e-12) ? Somin : 1e-6;                 // Small default to avoid division issues
      }

      // Step 2: Compute F_known from Mode 2 (imbibition scanning curve) at S_star
      // Mode 2 formula: F = (1 / (S - H + E) - 1 / E) / (1 / (Swma - H + E) - 1 / E)
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( wettingCurve,
                                                                  H,
                                                                  landParam,
                                                                  m_jerauldParam_a,
                                                                  m_jerauldParam_b,
                                                                  Scrt );
      real64 const Swma = 1 - (1 - Scrt);               // For two-phase, following existing Mode 2 pattern
      real64 F_known = 0.0;
      if( S_star <= Swma + 1e-12 )
      {
        real64 const F_num = (1. / (S_star - H + E) - 1. / E);
        real64 const F_denom = (1. / (Swma - H + E) - 1. / E);
        if( LvArray::math::abs( F_denom ) > 1e-12 )
        {
          F_known = F_num / F_denom;
        }
        else
        {
          F_known = (S_star > H) ? 1.0 : 0.0;
        }
      }
      else
      {
        // S_star > Swma: beyond Mode 2 range, use F_known = 0 (pure imbibition)
        F_known = 0.0;
      }
      F_known = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_known ));

      // Step 3: Compute F_star_target = 1 - F_known (from Eq. A-6)
      real64 const F_star_target = 1.0 - F_known;

      // Step 4: Solve quadratic for x = Sw_star_Hyst_ghost
      // F* = [1/(x - S_star + E) - 1/E] / [1/(x - H + E) - 1/E] = F_star_target
      // Note: Denominator uses H (first reversal point), not Sw_wr, so curve rejoins drainage at H
      // Rearranging gives: (1-T)*x^2 + B*x + C = 0
      // where T = F_star_target, a = S_star, b = H, e = E
      real64 const a = S_star;
      real64 const b = H;                // Use H instead of Sw_wr so curve rejoins drainage at H
      real64 const e = E;
      real64 const T = F_star_target;

      real64 const coeff_A = (1.0 - T);
      real64 const coeff_B = -(1.0 - T) * (a + b - 2.0 * e) - (1.0 - T) * e;
      real64 const coeff_C = (1.0 - T) * (a - e) * (b - e)
                             - e * e * (1.0 + T)
                             - e * (T * a - b);

      real64 x = S_star;
      real64 discriminant = coeff_B * coeff_B - 4.0 * coeff_A * coeff_C;

      if( discriminant >= 0.0 && LvArray::math::abs( coeff_A ) > 1e-14 )
      {
        real64 sqrt_disc = LvArray::math::sqrt( discriminant );
        real64 x1 = (-coeff_B + sqrt_disc) / (2.0 * coeff_A);
        real64 x2 = (-coeff_B - sqrt_disc) / (2.0 * coeff_A);

        // Pick the physically meaningful root: x must be > S_star
        if( x1 > S_star - 1e-8 && x2 > S_star - 1e-8 )
        {
          x = LvArray::math::min( x1, x2 );                 // Take the closer one
        }
        else if( x1 > S_star - 1e-8 )
        {
          x = x1;
        }
        else if( x2 > S_star - 1e-8 )
        {
          x = x2;
        }
      }

      // Step 5: Compute F_star for current saturation S
      // F_star = [1/(x - S + E) - 1/E] / [1/(x - H + E) - 1/E]
      // When S = H, F_star = 1.0 (pure drainage), so curve rejoins drainage at H
      real64 F_star = 1.0;
      real64 F_star_num = (1. / (x - S + E) - 1. / E);
      real64 F_star_denom = (1. / (x - H + E) - 1. / E);
      if( LvArray::math::abs( F_star_denom ) > 1e-12 )
      {
        F_star = F_star_num / F_star_denom;
        F_star = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_star ));
      }

      // Step 6: Compute Pc using Eq. A-5: P_c = P_c^Im + F* [P_c^Dr - P_c^Im]
      phaseCapPressure = pci + F_star * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F_star * (dpcd_dS - dpci_dS);
    }
    else
    {
      GEOS_THROW( GEOS_FMT( "{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                            "TableCapillaryPressureHysteresis",
                            (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                               ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                  InputError );
    }


  }

}
void TableCapillaryPressureHysteresis::KernelWrapper::computeTwoPhaseWetting( const geos::integer ipWetting,
                                                                              const geos::integer GEOS_UNUSED_PARAM( ipNonWetting ),
                                                                              const arraySlice1d< const geos::real64,
                                                                                                  compflow::USD_PHASE -
                                                                                                  1 > & phaseVolFraction,
                                                                              const arraySlice1d< const geos::real64,
                                                                                                  compflow::USD_PHASE
                                                                                                  -
                                                                                                  1 > & phaseMaxHistoricalVolFraction,
                                                                              const arraySlice1d< const geos::real64,
                                                                                                  compflow::USD_PHASE
                                                                                                  -
                                                                                                  1 > & phaseMinHistoricalVolFraction,
                                                                              const arraySlice1d< geos::real64,
                                                                                                  relperm::USD_RELPERM
                                                                                                  - 2 > & phaseTrappedVolFrac,
                                                                              arraySlice1d< geos::real64,
                                                                                            relperm::USD_RELPERM
                                                                                            -
                                                                                            2 > const & phaseCapPressure,
                                                                              arraySlice2d< geos::real64,
                                                                                            relperm::USD_RELPERM_DS
                                                                                            -
                                                                                            2 > const & dPhaseCapPressure_dPhaseVolFrac,
                                                                              ModeIndexType & mode,
                                                                              arraySlice1d< geos::real64,
                                                                                            compflow::USD_PHASE - 1 > & phaseMode2PeakVolFraction ) const
{
  using TTP = ThreePhasePairPhaseType;

  // Validate array sizes and indices before accessing
  GEOS_ASSERT_MSG( ipWetting >= 0, "ipWetting must be non-negative" );
  GEOS_ASSERT_MSG( static_cast< integer >(phaseVolFraction.size()) > ipWetting,
                   GEOS_FMT( "phaseVolFraction array too small: size={}, ipWetting={}. "
                             "This usually means the arrays haven't been properly resized. "
                             "Ensure resizeFields() has been called before using the KernelWrapper.",
                             phaseVolFraction.size(), ipWetting ));
  GEOS_ASSERT_MSG( static_cast< integer >(phaseMaxHistoricalVolFraction.size()) > ipWetting,
                   GEOS_FMT( "phaseMaxHistoricalVolFraction array too small: size={}, ipWetting={}. "
                             "This usually means the arrays haven't been properly resized. "
                             "Ensure resizeFields() has been called before using the KernelWrapper.",
                             phaseMaxHistoricalVolFraction.size(), ipWetting ));
  GEOS_ASSERT_MSG( static_cast< integer >(phaseMinHistoricalVolFraction.size()) > ipWetting,
                   GEOS_FMT( "phaseMinHistoricalVolFraction array too small: size={}, ipWetting={}. "
                             "This usually means the arrays haven't been properly resized. "
                             "Ensure resizeFields() has been called before using the KernelWrapper.",
                             phaseMinHistoricalVolFraction.size(), ipWetting ));
  GEOS_ASSERT_MSG( static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) >= 2,
                   GEOS_FMT( "m_wettingNonWettingCapillaryPressureKernelWrappers must have at least 2 elements, but got {}. "
                             "This usually means createAllTableKernelWrappers() failed to populate the arrays.",
                             m_wettingNonWettingCapillaryPressureKernelWrappers.size()));

  // TEMPORARY: Disable scanning curves for testing convergence
  // Define at function scope so it's accessible to both wetting and non-wetting phase code
  constexpr bool ENABLE_SCANNING_CURVES = true;

  // TEMPORARY: Force IMBIBITION mode for testing (overrides normal mode determination)
  constexpr bool FORCE_IMBIBITION_MODE = false;

  // TEMPORARY: Enable/disable Mode 4 (IMBIBITION_TO_DRAINAGE_FROM_SCANNING)
  // When false, the code stays in Mode 2 and does not transition to Mode 4
  constexpr bool ENABLE_MODE4 = false;

  // Determine mode based on saturation condition and flow direction
  // Use DRAINAGE when saturation is at or below minimum historical
  // Use scanning curves when above minimum, detecting direction of change
  bool const useDrainage = FORCE_IMBIBITION_MODE ? false :
                           (!m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] ||
                            phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer);

  //--- wetting  cap pressure -- W/O or W/G two phase flow
  // Use drainage curve when S <= S_min
  // Use scanning curves when S > S_min, but detect if we're in secondary drainage
  // DEBUG: Print mode for capillary pressure
  // printf("CapPressure: mode=%d, S_w=%.6e, S_min=%.6e, S_max=%.6e, hasHyst=%d\n",
  //        static_cast<int>(mode),
  //        phaseVolFraction[ipWetting],
  //        phaseMinHistoricalVolFraction[ipWetting],
  //        phaseMaxHistoricalVolFraction[ipWetting],
  //        static_cast<int>(m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING]));

  if( useDrainage )
  {
    // Use simple drainage curve (matching relative permeability)
    mode = ModeIndexType::DRAINAGE;
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting],
                                                         m_wettingCurve.oppositeBoundPhaseVolFraction );
    // printf("CapPressure: Using DRAINAGE curve (arrayIndex=0)\n");
    GEOS_ASSERT_MSG( static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) > ModeIndexType::DRAINAGE,
                     "Invalid array index for kernel wrapper access" );
    computeBoundCapillaryPressure(
      m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE],
      phaseVolFraction[ipWetting],
      phaseCapPressure[ipWetting],
      dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting] );
  }
  else
  {
    // Use scanning curves - original conservative approach: only switch from pure states to scanning curves
    // Scanning curve modes persist once set (no switching between DRAINAGE_TO_IMBIBITION and IMBIBITION_TO_DRAINAGE)

    // If mode is already set to pure IMBIBITION, use it directly (for Pc bounds computation)
    if( mode == ModeIndexType::IMBIBITION )
    {
      // Force pure imbibition mode - use the imbibition table directly
      phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting],
                                                           m_wettingCurve.oppositeBoundPhaseVolFraction );
      GEOS_ASSERT_MSG( static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) > ModeIndexType::IMBIBITION,
                       "Invalid array index for kernel wrapper access" );
      computeBoundCapillaryPressure(
        m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::IMBIBITION],
        phaseVolFraction[ipWetting],
        phaseCapPressure[ipWetting],
        dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting] );
      return;               // Exit early, don't use scanning curves
    }

    // If FORCE_IMBIBITION_MODE is enabled, override mode to IMBIBITION
    if( FORCE_IMBIBITION_MODE )
    {
      mode = ModeIndexType::IMBIBITION;
    }

    // If scanning curves are disabled, reset any scanning curve modes back to pure modes
    if( !ENABLE_SCANNING_CURVES )
    {
      if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION ||
          mode == ModeIndexType::IMBIBITION_TO_DRAINAGE ||
          mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
      {
        // Reset to pure drainage or imbibition based on saturation
        if( FORCE_IMBIBITION_MODE )
        {
          mode = ModeIndexType::IMBIBITION;
        }
        else
        {
          real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
          real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
          real64 const Shy_min = (phaseMinHistoricalVolFraction[ipWetting] > Smin) ?
                                 phaseMinHistoricalVolFraction[ipWetting] : Smin;
          real64 const Shy_max = (phaseMaxHistoricalVolFraction[ipWetting] < Smax) ?
                                 phaseMaxHistoricalVolFraction[ipWetting] : Smax;
          real64 const currentS = phaseVolFraction[ipWetting];

          if( currentS <= Shy_min + flowReversalBuffer )
          {
            mode = ModeIndexType::DRAINAGE;
          }
          else if( currentS >= Shy_max - flowReversalBuffer )
          {
            mode = ModeIndexType::IMBIBITION;
          }
          else
          {
            // Default to drainage if in between
            mode = ModeIndexType::DRAINAGE;
          }
        }
      }
    }

    bool justSwitchedToMode2 = false;
    if( ENABLE_SCANNING_CURVES && (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION))
    {
      // Transitioning from a pure state to a scanning curve.
      // The previous mode tells us the direction of the reversal:
      //   - From DRAINAGE: saturation is now increasing → DRAINAGE_TO_IMBIBITION (Mode 2)
      //   - From IMBIBITION: saturation is now decreasing → IMBIBITION_TO_DRAINAGE (Mode 3)
      if( mode == ModeIndexType::DRAINAGE )
      {
        mode = ModeIndexType::DRAINAGE_TO_IMBIBITION;
        justSwitchedToMode2 = true;
      }
      else
      {
        mode = ModeIndexType::IMBIBITION_TO_DRAINAGE;
        // Reset Mode 2 peak when leaving imbibition
        if( phaseMode2PeakVolFraction[ipWetting] > 0.0 )
        {
          phaseMode2PeakVolFraction[ipWetting] = 0.0;
        }
      }
    }
    // Mode 2 peak is now tracked explicitly at end of timestep in saveConvergedPhaseVolFractionState()
    // No need to update it during Newton iterations - use constant value from previous timestep
    // This prevents oscillations during Newton solve from affecting the peak value
    // Check for transition from DRAINAGE_TO_IMBIBITION to drainage (departing from scanning curve)
    // Mode 4: Switch when saturation decreases from the peak reached during Mode 2
    // In Mode 2, saturation increases from Shy_min (reversal point) towards Mode2_peak
    // Once saturation starts decreasing from Mode2_peak by at least 1e-4, switch to Mode 4
    // IMPORTANT: Only check for Mode 4 if we were already in Mode 2 (not just switched from Mode 0)
    // This prevents immediate Mode 0 -> Mode 2 -> Mode 4 transition within a single Newton iteration
    // TEMPORARY: Disable scanning curves for testing convergence
    if( ENABLE_MODE4 && ENABLE_SCANNING_CURVES && mode == ModeIndexType::DRAINAGE_TO_IMBIBITION && !justSwitchedToMode2 )
    {
      real64 const currentS = phaseVolFraction[ipWetting];
      real64 const Mode2_peak = phaseMode2PeakVolFraction[ipWetting];
      // Use same approach as other mode switches: small buffer + mode persistence
      // Similar to Mode 0 -> Mode 2: use flowReversalBuffer (1e-12) for consistency
      // Mode 4 will persist once set, preventing oscillation-induced switching
      // However, we need a slightly larger buffer than 1e-12 because:
      // 1. The peak is a dynamic value that was just reached
      // 2. Saturation can oscillate around the peak during Newton iterations
      // Use a relative buffer: 0.01% of the peak value, with a minimum of 1e-4
      // This is still much smaller than before but provides some protection against numerical noise
      real64 const relativeBuffer = 0.0001 * Mode2_peak;               // 0.01% of peak
      real64 const absoluteBuffer = 1e-4;               // Minimum absolute buffer (same as original)
      real64 const decreaseBuffer = LvArray::math::max( relativeBuffer, absoluteBuffer );
      // Valid range for Mode 4 switching: between min and max historical saturations
      real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
      real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
      real64 const S_test_min = (phaseMinHistoricalVolFraction[ipWetting] > Smin) ?
                                phaseMinHistoricalVolFraction[ipWetting] : Smin;
      real64 const S_test_max = (phaseMaxHistoricalVolFraction[ipWetting] < Smax) ?
                                phaseMaxHistoricalVolFraction[ipWetting] : Smax;

      // Check: saturation must have decreased from the Mode 2 peak
      // Only switch to Mode 4 if currentS < Mode2_peak - decreaseBuffer
      // Also ensure Mode2_peak is valid (non-zero, meaning we've tracked a peak)
      // And current saturation is within valid historical range
      // Once Mode 4 is set, it persists (like Mode 2), preventing oscillation-induced switching
      bool saturationHasDecreased = (Mode2_peak > 0.0 && currentS < Mode2_peak - decreaseBuffer);

      if( saturationHasDecreased &&
          currentS >= S_test_min &&
          currentS <= S_test_max )
      {
        mode = ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING;
      }
      // Otherwise, keep Mode 2 (saturation still increasing or hasn't decreased enough from peak)
      // Mode 2 persists once set, similar to how Mode 0 persists when S <= S_min + flowReversalBuffer
    }
    // IMBIBITION_TO_DRAINAGE_FROM_SCANNING mode persists once set

    // If scanning curves are disabled, use pure drainage/imbibition curves
    if( !ENABLE_SCANNING_CURVES &&
        (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION))
    {
      // Use pure drainage or imbibition curve
      phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting],
                                                           m_wettingCurve.oppositeBoundPhaseVolFraction );
      GEOS_ASSERT_MSG( static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()) > mode,
                       "Invalid array index for kernel wrapper access" );
      computeBoundCapillaryPressure(
        m_wettingNonWettingCapillaryPressureKernelWrappers[mode],
        phaseVolFraction[ipWetting],
        phaseCapPressure[ipWetting],
        dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting] );
    }
    else
    {
      // Use scanning curves
      // printf("CapPressure: Using scanning curve (mode=%d)\n", static_cast<int>(mode));
      computeImbibitionWettingCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                 m_wettingCurve,
                                                 m_landParam[ipWetting],
                                                 phaseVolFraction[ipWetting],
                                                 phaseMinHistoricalVolFraction[ipWetting],
                                                 phaseMaxHistoricalVolFraction[ipWetting],
                                                 phaseMode2PeakVolFraction[ipWetting],
                                                 phaseTrappedVolFrac[ipWetting],
                                                 phaseCapPressure[ipWetting],
                                                 dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting],
                                                 mode );
    }
  }

// trapped vol fraction
  if( mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION ||
      mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
  {



    real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
    real64 const Shy = (phaseMinHistoricalVolFraction[ipWetting] < Smin)
                                       ? phaseMinHistoricalVolFraction[ipWetting]
                                       : Smin;
    real64 Scrt = 0.0;
    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, Shy,
                                                                m_landParam[ipWetting],
                                                                m_jerauldParam_a,
                                                                m_jerauldParam_b,
                                                                Scrt );
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );


    //keep the same Land coeff as two phase only
    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve.toNonWetting(), Shy,
                                                                m_landParam[ipWetting],
                                                                m_jerauldParam_a,
                                                                m_jerauldParam_b,
                                                                Scrt );
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );



  }
  else if( mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
  {

    real64 const Smax = m_wettingCurve.imbibitionExtremaPhaseVolFraction;
    real64 const Shy = (phaseMaxHistoricalVolFraction[ipWetting] < Smax)
                                       ? phaseMaxHistoricalVolFraction[ipWetting]
                                       : Smax;
    real64 Scrt = 0.0;
    //TODO (jacques) check if still accurate
    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve,
                                                                Shy,
                                                                m_landParam[ipWetting],
                                                                m_jerauldParam_a,
                                                                m_jerauldParam_b,
                                                                Scrt );

    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );

    //keep the same Land coeff as two phase only
    KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve.toNonWetting(), Shy,
                                                                m_landParam[ipWetting],
                                                                m_jerauldParam_a,
                                                                m_jerauldParam_b,
                                                                Scrt );
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );

  }

}

void TableCapillaryPressureHysteresis::KernelWrapper::computeTwoPhaseNonWetting( const geos::integer ipWetting,
                                                                                 const geos::integer ipNonWetting,
                                                                                 const arraySlice1d< const geos::real64,
                                                                                                     compflow::USD_PHASE -
                                                                                                     1 > & phaseVolFraction,
                                                                                 const arraySlice1d< const geos::real64,
                                                                                                     compflow::USD_PHASE
                                                                                                     -
                                                                                                     1 > & phaseMaxHistoricalVolFraction,
                                                                                 const arraySlice1d< const geos::real64,
                                                                                                     compflow::USD_PHASE
                                                                                                     -
                                                                                                     1 > & phaseMinHistoricalVolFraction,
                                                                                 const arraySlice1d< geos::real64,
                                                                                                     relperm::USD_RELPERM
                                                                                                     -
                                                                                                     2 > & phaseTrappedVolFrac,
                                                                                 arraySlice1d< geos::real64,
                                                                                               relperm::USD_RELPERM
                                                                                               -
                                                                                               2 > const & phaseCapPressure,
                                                                                 arraySlice2d< geos::real64,
                                                                                               relperm::USD_RELPERM_DS
                                                                                               -
                                                                                               2 > const & dPhaseCapPressure_dPhaseVolFrac,
                                                                                 ModeIndexType & mode ) const
{
  using TTP = ThreePhasePairPhaseType;

  // TEMPORARY: Disable scanning curves for testing convergence
  // Define at function scope so it's accessible throughout this function
  constexpr bool ENABLE_SCANNING_CURVES = true;

  // TEMPORARY: Force IMBIBITION mode for testing (overrides normal mode determination)
  constexpr bool FORCE_IMBIBITION_MODE = false;

  //update state
  // TODO check if we can get rid of  DRAINAGE_TO_IMBIBITION && IMBIBITION_TO_DRAINAGE
  if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION &&
      phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
    mode = ModeIndexType::DRAINAGE;
  if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE &&
      phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
    mode = ModeIndexType::IMBIBITION;



  // If scanning curves are disabled, reset any scanning curve modes back to pure modes
  if( !ENABLE_SCANNING_CURVES )
  {
    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION ||
        mode == ModeIndexType::IMBIBITION_TO_DRAINAGE ||
        mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
    {
      // Reset to pure drainage or imbibition based on saturation
      if( FORCE_IMBIBITION_MODE )
      {
        mode = ModeIndexType::IMBIBITION;
      }
      else
      {
        if( phaseVolFraction[ipNonWetting] <= phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
        {
          mode = ModeIndexType::DRAINAGE;
        }
        else if( phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] - flowReversalBuffer )
        {
          mode = ModeIndexType::IMBIBITION;
        }
        else
        {
          // Default to drainage if in between
          mode = ModeIndexType::DRAINAGE;
        }
      }
    }
  }

  // Force IMBIBITION mode if flag is set (overrides all other logic)
  if( FORCE_IMBIBITION_MODE )
  {
    mode = ModeIndexType::IMBIBITION;
  }

  if( ENABLE_SCANNING_CURVES && (mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION))
  {
    // Handle transitions from pure states to scanning curves
    // If current saturation is significantly above min historical, we're in imbibition
    if( phaseVolFraction[ipNonWetting] > phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
    {
      if( mode == ModeIndexType::DRAINAGE )
      {
        mode = ModeIndexType::DRAINAGE_TO_IMBIBITION;
      }
    }
    // If current saturation is significantly below max historical, we're in drainage
    else if( phaseVolFraction[ipNonWetting] < phaseMaxHistoricalVolFraction[ipNonWetting] - flowReversalBuffer )
    {
      if( mode == ModeIndexType::IMBIBITION )
      {
        mode = ModeIndexType::IMBIBITION_TO_DRAINAGE;
      }
    }
  }
  // Mode 4 switching is only implemented for wetting phase, not for non-wetting phase
  // Scanning curve modes persist - no switching between DRAINAGE_TO_IMBIBITION and IMBIBITION_TO_DRAINAGE

  // Use simple drainage/imbibition curves when:
  // 1. No hysteresis enabled, OR
  // 2. We're in pure DRAINAGE mode (use drainage curve), OR
  // 3. We're in pure IMBIBITION mode (use imbibition curve)
  // Use scanning curves only during transitions (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
  if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] ||
      mode == ModeIndexType::DRAINAGE ||
      mode == ModeIndexType::IMBIBITION )
  {
    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting],
                                                            (mode == ModeIndexType::DRAINAGE)
                                                                       ? m_nonWettingCurve.drainageExtremaPhaseVolFraction
                                                                       : m_nonWettingCurve.imbibitionExtremaPhaseVolFraction );
    // Ensure mode is a valid array index (0 or 1)
    integer const arrayIndex = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE : ModeIndexType::IMBIBITION;
    GEOS_ASSERT_MSG( arrayIndex >= 0 && arrayIndex < static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()),
                     "Invalid array index for kernel wrapper access" );
    computeBoundCapillaryPressure(
      m_wettingNonWettingCapillaryPressureKernelWrappers[arrayIndex],
      phaseVolFraction[ipNonWetting],
      phaseCapPressure[ipNonWetting],
      dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *= -1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;

  }
  else
  {
    // We're in a transition state (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
    // Use scanning curves (unless disabled)
    if( !ENABLE_SCANNING_CURVES )
    {
      // Scanning curves disabled - use pure drainage/imbibition based on current mode
      // This should not happen if reset logic worked correctly, but handle it anyway
      if( mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION )
      {
        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting],
                                                                (mode == ModeIndexType::DRAINAGE)
                                                                              ? m_nonWettingCurve.drainageExtremaPhaseVolFraction
                                                                              : m_nonWettingCurve.imbibitionExtremaPhaseVolFraction );
        integer const arrayIndex = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE : ModeIndexType::IMBIBITION;
        GEOS_ASSERT_MSG( arrayIndex >= 0 && arrayIndex < static_cast< integer >(m_wettingNonWettingCapillaryPressureKernelWrappers.size()),
                         "Invalid array index for kernel wrapper access" );
        computeBoundCapillaryPressure(
          m_wettingNonWettingCapillaryPressureKernelWrappers[arrayIndex],
          phaseVolFraction[ipNonWetting],
          phaseCapPressure[ipNonWetting],
          dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
        phaseCapPressure[ipNonWetting] *= -1;
        dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
      }
      else
      {
        mode = ModeIndexType::DRAINAGE;
        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting],
                                                                m_nonWettingCurve.drainageExtremaPhaseVolFraction );
        computeBoundCapillaryPressure(
          m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE],
          phaseVolFraction[ipNonWetting],
          phaseCapPressure[ipNonWetting],
          dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
        phaseCapPressure[ipNonWetting] *= -1;
        dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
      }
    }
    else
    {
      // Use scanning curves

      computeImbibitionNonWettingCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                    m_nonWettingCurve,
                                                    m_landParam[ipNonWetting],
                                                    phaseVolFraction[ipNonWetting],
                                                    phaseMaxHistoricalVolFraction[ipNonWetting],
                                                    phaseTrappedVolFrac[ipNonWetting],
                                                    phaseCapPressure[ipNonWetting],
                                                    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting],
                                                    mode );

      // when pc is on the gas phase, we need to multiply user input by -1
      // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
      phaseCapPressure[ipNonWetting] *= -1;
      dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
    }
  }

// trapped vol fraction
  if( mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION ||
      mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
  {

    {
      real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
      real64 const Shy = (phaseMaxHistoricalVolFraction[ipNonWetting] > Smax)
                                       ? phaseMaxHistoricalVolFraction[ipNonWetting]
                                       : Smax;
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, Shy,
                                                                  m_landParam[ipNonWetting],
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );

      //keep the same Land coeff as two phase only
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve.toWetting(), Shy,
                                                                  m_landParam[ipNonWetting],
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );
    }
  }
  else if( mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
  {

    {
      real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;;
      real64 const Shy = (phaseMinHistoricalVolFraction[ipNonWetting] > Smin)
                                       ? phaseMinHistoricalVolFraction[ipNonWetting]
                                       : Smin;
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, Shy,
                                                                  m_landParam[ipNonWetting],
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );

      //keep the same Land coeff as two phase only
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve.toWetting(), Shy,
                                                                  m_landParam[ipNonWetting],
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );
    }
  }


}

void TableCapillaryPressureHysteresis::KernelWrapper::computeThreePhase( const geos::integer ipWetting,
                                                                         const geos::integer GEOS_UNUSED_PARAM( ipInter ),
                                                                         const geos::integer ipNonWetting,
                                                                         const arraySlice1d< const geos::real64,
                                                                                             compflow::USD_PHASE -
                                                                                             1 > & phaseVolFraction,
                                                                         const arraySlice1d< const geos::real64,
                                                                                             compflow::USD_PHASE
                                                                                             -
                                                                                             1 > & phaseMaxHistoricalVolFraction,
                                                                         const arraySlice1d< const geos::real64,
                                                                                             compflow::USD_PHASE
                                                                                             -
                                                                                             1 > & phaseMinHistoricalVolFraction,
                                                                         const arraySlice1d< geos::real64,
                                                                                             relperm::USD_RELPERM
                                                                                             - 2 > & phaseTrappedVolFrac,
                                                                         const arraySlice1d< geos::real64,
                                                                                             relperm::USD_RELPERM -
                                                                                             2 > & phaseCapPressure,
                                                                         const arraySlice2d< geos::real64,
                                                                                             relperm::USD_RELPERM_DS
                                                                                             -
                                                                                             2 > & dPhaseCapPressure_dPhaseVolFrac,
                                                                         ModeIndexType & mode,
                                                                         arraySlice1d< geos::real64,
                                                                                       compflow::USD_PHASE
                                                                                       -
                                                                                       1 > & phaseMode2PeakVolFraction ) const
{


  LvArray::forValuesInSlice( dPhaseCapPressure_dPhaseVolFrac, []( real64 & val ) { val = 0.0; } );
  using TTP = ThreePhasePairPhaseType;

  // -- wetting curve if drainage only
  constexpr bool ENABLE_COMPUTE_BRANCH_DEBUG = false;
  bool usePureCurve = !m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] ||
                      (mode == ModeIndexType::DRAINAGE &&
                       phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer) ||
                      (mode == ModeIndexType::IMBIBITION &&
                       phaseVolFraction[ipWetting] >= phaseMaxHistoricalVolFraction[ipWetting] + flowReversalBuffer);

  if constexpr (ENABLE_COMPUTE_BRANCH_DEBUG) {
    if( mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::IMBIBITION )
    {
      std::cout << "[COMPUTE_BRANCH_DEBUG] Two-phase wetting, mode=" << (mode == ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                << ", S=" << phaseVolFraction[ipWetting]
                << ", S_min=" << phaseMinHistoricalVolFraction[ipWetting]
                << ", S_max=" << phaseMaxHistoricalVolFraction[ipWetting]
                << ", usePureCurve=" << (usePureCurve ? "true" : "false")
                << ", hasHysteresis=" << m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] << std::endl;
    }
  }

  if( usePureCurve )
  {
    // water-oil capillary pressure
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting],
                                                         m_wettingCurve.oppositeBoundPhaseVolFraction );
    constexpr bool ENABLE_COMPUTE_TABLE_DEBUG = false;
    real64 S_before = phaseVolFraction[ipWetting];
    real64 dPc_dS_before = dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting];
    integer const tableIndex = static_cast< integer >(mode);
    if constexpr (ENABLE_COMPUTE_TABLE_DEBUG) {
      std::cout << "[COMPUTE_TABLE_DEBUG_TCH] Forward table lookup: mode=" << (mode == ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                << ", tableIndex=" << tableIndex
                << ", S=" << S_before
                << ", wrapper.size()=" << m_wettingIntermediateCapillaryPressureKernelWrappers.size() << std::endl;
    }
    phaseCapPressure[ipWetting] =
      m_wettingIntermediateCapillaryPressureKernelWrappers[mode].compute(
        &(phaseVolFraction)[ipWetting],
        &(dPhaseCapPressure_dPhaseVolFrac)[ipWetting][ipWetting] );
    real64 dPc_dS_after = dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting];
    if constexpr (ENABLE_COMPUTE_TABLE_DEBUG) {
      std::cout << "[COMPUTE_TABLE_DEBUG_TCH] After table.compute: Pc=" << phaseCapPressure[ipWetting]
                << ", dPc_dS_before=" << dPc_dS_before
                << ", dPc_dS_after=" << dPc_dS_after << std::endl;
    }
  }
  else
  {
    // We're in a scanning curve state - preserve the mode if already set, otherwise determine from position
    if( mode != ModeIndexType::DRAINAGE_TO_IMBIBITION &&
        mode != ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      // Starting from a pure state - determine initial scanning curve direction
      mode = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE_TO_IMBIBITION
                                                         : ModeIndexType::IMBIBITION_TO_DRAINAGE;
    }
    computeImbibitionWettingCapillaryPressure( m_wettingIntermediateCapillaryPressureKernelWrappers,
                                               m_wettingCurve,
                                               m_nonWettingCurve,
                                               m_landParam[ipWetting],
                                               phaseVolFraction[ipWetting],
                                               phaseMinHistoricalVolFraction[ipWetting],
                                               phaseMaxHistoricalVolFraction[ipWetting],
                                               phaseMode2PeakVolFraction[ipWetting],
                                               phaseTrappedVolFrac[ipWetting],
                                               phaseCapPressure[ipWetting],
                                               dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting],
                                               mode );


  }


  // -- non-wetting cure if drainage only
  // gas-oil capillary pressure
  if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] ||
      (mode == ModeIndexType::DRAINAGE &&
       phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer) ||
      (mode == ModeIndexType::IMBIBITION &&
       phaseVolFraction[ipNonWetting] <= phaseMinHistoricalVolFraction[ipNonWetting] + flowReversalBuffer))
  {
    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting],
                                                            (mode == ModeIndexType::DRAINAGE)
                                                                       ? m_nonWettingCurve.drainageExtremaPhaseVolFraction
                                                                       : m_nonWettingCurve.imbibitionExtremaPhaseVolFraction );
    phaseCapPressure[ipNonWetting] =
      m_nonWettingIntermediateCapillaryPressureKernelWrappers[mode].compute(
        &(phaseVolFraction)[ipNonWetting],
        &(dPhaseCapPressure_dPhaseVolFrac)[ipNonWetting][ipNonWetting] );


    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *= -1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
  }
  else
  {
    // We're in a scanning curve state - preserve the mode if already set, otherwise determine from position
    if( mode != ModeIndexType::DRAINAGE_TO_IMBIBITION &&
        mode != ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      // Starting from a pure state - determine initial scanning curve direction
      mode = (mode == ModeIndexType::DRAINAGE) ? ModeIndexType::DRAINAGE_TO_IMBIBITION
                                                         : ModeIndexType::IMBIBITION_TO_DRAINAGE;
    }

    computeImbibitionNonWettingCapillaryPressure( m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                                                  m_nonWettingCurve,
                                                  m_wettingCurve,
                                                  m_landParam[ipNonWetting],
                                                  phaseVolFraction[ipNonWetting],
                                                  phaseMinHistoricalVolFraction[ipNonWetting],
                                                  phaseTrappedVolFrac[ipNonWetting],
                                                  phaseCapPressure[ipNonWetting],
                                                  dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting],
                                                  mode );
    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *= -1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;

    //update trapped fraction
    if( mode == ModeIndexType::DRAINAGE || mode == ModeIndexType::DRAINAGE_TO_IMBIBITION ||
        mode == ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
    {


      {
        real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
        real64 const Shy = (phaseMinHistoricalVolFraction[ipWetting] < Smin)
                                           ? phaseMinHistoricalVolFraction[ipWetting]
                                           : Smin;
        real64 Scrt = 0.0;
        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, Shy,
                                                                    m_landParam[ipWetting],
                                                                    m_jerauldParam_a, m_jerauldParam_b,
                                                                    Scrt );
        phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );
      }

      {
        real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
        real64 const Shy = (phaseMaxHistoricalVolFraction[ipNonWetting] > Smax)
                                           ? phaseMaxHistoricalVolFraction[ipNonWetting]
                                           : Smax;
        real64 Scrt = 0.0;
        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, Shy,
                                                                    m_landParam[ipNonWetting],
                                                                    m_jerauldParam_a, m_jerauldParam_b,
                                                                    Scrt );
        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );
      }
    }
    else if( mode == ModeIndexType::IMBIBITION || mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      {
        real64 const Smax = m_wettingCurve.imbibitionExtremaPhaseVolFraction;
        real64 const Shy = (phaseMaxHistoricalVolFraction[ipWetting] < Smax)
                                           ? phaseMaxHistoricalVolFraction[ipWetting]
                                           : Smax;
        real64 Scrt = 0.0;
        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, Shy,
                                                                    m_landParam[ipWetting],
                                                                    m_jerauldParam_a, m_jerauldParam_b,
                                                                    Scrt );
        phaseTrappedVolFrac[ipWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipWetting] );
      }

      {
        real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;;
        real64 const Shy = (phaseMinHistoricalVolFraction[ipNonWetting] > Smin)
                                           ? phaseMinHistoricalVolFraction[ipNonWetting]
                                           : Smin;
        real64 Scrt = 0.0;
        KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, Shy,
                                                                    m_landParam[ipNonWetting],
                                                                    m_jerauldParam_a, m_jerauldParam_b,
                                                                    Scrt );
        phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( Scrt, phaseVolFraction[ipNonWetting] );
      }
    }


  }


}


void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionNonWettingCapillaryPressure(
  const arrayView1d< const TableFunction::KernelWrapper > & nonWettingKernelWrapper,
  const KilloughHysteresis::HysteresisCurve & nonWettingCurve,
  const KilloughHysteresis::HysteresisCurve & wettingCurve,
  const geos::real64 & landParam,
  const geos::real64 & phaseVolFraction,
  const geos::real64 & phaseMaxHistoricalVolFraction,
  geos::real64 & phaseTrappedVolFrac,
  geos::real64 & phaseCapPressure,
  geos::real64 & dPhaseCapPressure_dPhaseVolFrac,
  const ModeIndexType & mode ) const
{
  GEOS_ASSERT( !nonWettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smii = nonWettingCurve.imbibitionExtremaPhaseVolFraction;
  real64 const Smid = nonWettingCurve.drainageExtremaPhaseVolFraction;
  real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;

  GEOS_UNUSED_VAR( Smii, Smid );

//  if( S >= Smax )
//  {
//    //above accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else if( S <= Smid )
//  {
//    //below accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = nonWettingKernelWrapper[ModeIndexType::IMBIBITION].compute( &S, &dpci_dS );
    real64 const pcd = nonWettingKernelWrapper[ModeIndexType::DRAINAGE].compute( &S, &dpcd_dS );

    // Step 1: get the trapped from wetting data
    real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

    //drainage to imbibition
    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
    {
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( nonWettingCurve, Shy, landParam,
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );

      real64 const E = m_killoughCurvatureParamCapPres;

      //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
      real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. compute dF_dS
      real64 dF_dS = (1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

      //Step 4. Eventually assemble everything following (EQ. 34.20)
      phaseCapPressure = pcd + F * (pci - pcd);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpci_dS - dpcd_dS);
      dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pci - pcd);

      //update trapped fraction
      phaseTrappedVolFrac = LvArray::math::min( Scrt, S );

    }
    //imbibition to drainage
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( wettingCurve, Shy, landParam,
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      real64 Sgma = 1. - Scrt - m_phaseIntermediateMinVolFraction;

      real64 const E = m_killoughCurvatureParamCapPres;

      //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
      real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Sgma - Shy + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. compute dF_dS
      real64 dF_dS = (-1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

      //Step 4. Eventually assemble everything following (EQ. 34.20)
      phaseCapPressure = pci + F * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpcd_dS - dpci_dS);
      dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pcd - pci);
    }
    else
    {
      GEOS_THROW( GEOS_FMT( "{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                            "TableCapillaryPressureHysteresis",
                            (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                               ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                  InputError );
    }


  }
}


void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionNonWettingCapillaryPressure(
  const arrayView1d< const TableFunction::KernelWrapper > & nonWettingKernelWrapper,
  const KilloughHysteresis::HysteresisCurve & nonWettingCurve,
  const geos::real64 & landParam,
  const geos::real64 & phaseVolFraction,
  const geos::real64 & phaseMaxHistoricalVolFraction,
  geos::real64 & phaseTrappedVolFrac,
  geos::real64 & phaseCapPressure,
  geos::real64 & dPhaseCapPressure_dPhaseVolFrac,
  const ModeIndexType & mode ) const
{

  GEOS_ASSERT( !nonWettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smii = nonWettingCurve.imbibitionExtremaPhaseVolFraction;
  real64 const Smid = nonWettingCurve.drainageExtremaPhaseVolFraction;
  real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;

  GEOS_UNUSED_VAR( Smii, Smid );

//  if( S >= Smax )
//  {
//    //above accessible range
//    phaseCapPressure = CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else if( S <= Smid )
//  {
//    //below accessible range
//    phaseCapPressure = -CAP_INF;
//    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
//  }
//  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = nonWettingKernelWrapper[ModeIndexType::IMBIBITION].compute( &S, &dpci_dS );
    real64 const pcd = nonWettingKernelWrapper[ModeIndexType::DRAINAGE].compute( &S, &dpcd_dS );

    // Step 1: get the trapped from wetting data
    real64 const Shy = (phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax;

    //drainage to imbibition
    if( mode == ModeIndexType::DRAINAGE_TO_IMBIBITION )
    {
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( nonWettingCurve, Shy, landParam,
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );

      real64 const E = m_killoughCurvatureParamCapPres;

      //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
      real64 F = (1. / (Shy - S + E) - 1. / E) / (1. / (Shy - Scrt + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. compute dF_dS
      real64 dF_dS = (1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

      //Step 4. Eventually assemble everything following (EQ. 34.20)
      phaseCapPressure = pcd + F * (pci - pcd);
      dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * (dpci_dS - dpcd_dS);
      dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pci - pcd);

      //update trapped fraction
      phaseTrappedVolFrac = LvArray::math::min( Scrt, S );

    }
    //imbibition to drainage
    else if( mode == ModeIndexType::IMBIBITION_TO_DRAINAGE )
    {
      real64 Scrt = 0.0;
      KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( nonWettingCurve, Shy, landParam,
                                                                  m_jerauldParam_a, m_jerauldParam_b,
                                                                  Scrt );
      real64 Sgma = 1. - (1. - Scrt);

      real64 const E = m_killoughCurvatureParamCapPres;

      //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
      real64 F = (1. / (S - Shy + E) - 1. / E) / (1. / (Sgma - Shy + E) - 1. / E);
      //force bound
      F = LvArray::math::max( F, 0.0 );
      F = LvArray::math::min( F, 1.0 );

      //Step 3. compute dF_dS
      real64 dF_dS = (-1. / (S * S)) / (1. / (Shy - Scrt + E) - 1. / E);

      //Step 4. Eventually assemble everything following (EQ. 34.20)
      phaseCapPressure = pci + F * (pcd - pci);
      dPhaseCapPressure_dPhaseVolFrac = dpcd_dS + F * (dpcd_dS - dpci_dS);
      dPhaseCapPressure_dPhaseVolFrac += dF_dS * (pcd - pci);
    }
    else
    {
      GEOS_THROW( GEOS_FMT( "{}: State is {}.Shouldnt be used in pure DRAINAGE or IMBIBITION.",
                            "TableCapillaryPressureHysteresis",
                            (mode == ModeIndexType::DRAINAGE) ? "DRAINAGE" : ((mode ==
                                                                               ModeIndexType::IMBIBITION)
                                                                                            ? "IMBIBITION"
                                                                                            : "UNKNOWN")),
                  InputError );
    }


  }
}


void TableCapillaryPressureHysteresis::KernelWrapper::computeBoundCapillaryPressure(
  const TableFunction::KernelWrapper & drainageRelpermWrapper,
  const geos::real64 & phaseVolFraction,
  geos::real64 & phaseCapPressure,
  geos::real64 & dPhaseCapPressure_dPhaseVolFrac ) const
{
  phaseCapPressure = drainageRelpermWrapper.compute( &phaseVolFraction,
                                                     &dPhaseCapPressure_dPhaseVolFrac );
}

/// Helper function to compute Pc(S) for given S, mode, and historical values
/// Used for Newton-Raphson inversion in computeInv
GEOS_HOST_DEVICE
inline real64 TableCapillaryPressureHysteresis::KernelWrapper::computeCapillaryPressureForSaturation(
  real64 const S,
  fields::cappres::ModeIndexType const & mode,
  integer const ipPhase,
  real64 const & phaseMinHistoricalVolFraction,
  real64 const & phaseMaxHistoricalVolFraction,
  real64 const & phaseMode2PeakVolFraction,
  arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
  KilloughHysteresis::HysteresisCurve const & wettingCurve,
  KilloughHysteresis::HysteresisCurve const & nonWettingCurve,
  real64 const & landParam,
  real64 const & phaseIntermediateMinVolFraction,
  real64 const & killoughCurvatureParam,
  real64 const & jerauldParam_a,
  real64 const & jerauldParam_b,
  bool const isWettingPhase,
  real64 const precomputedScrt,
  real64 const precomputedDenomF,
  real64 const precomputedShy ) const
{
  GEOS_UNUSED_VAR( ipPhase, landParam, jerauldParam_a, jerauldParam_b );
  real64 pc = 0.0;
  real64 dpc_dS = 0.0;

  // For pure drainage or imbibition modes, use the table directly
  if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
  {
    integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
    if( arrayIndex < static_cast< integer >(capPresKernelWrappers.size()))
    {
      constexpr bool ENABLE_COMPUTE_PC_DEBUG = false;
      if constexpr (ENABLE_COMPUTE_PC_DEBUG) {
        std::cout << "[COMPUTE_PC_DEBUG] computeCapillaryPressureForSaturation: mode="
                  << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                  << ", arrayIndex=" << arrayIndex
                  << ", S=" << S
                  << ", capPresKernelWrappers.size()=" << capPresKernelWrappers.size() << std::endl;
      }
      pc = capPresKernelWrappers[arrayIndex].compute( &S, &dpc_dS );
      if constexpr (ENABLE_COMPUTE_PC_DEBUG) {
        std::cout << "[COMPUTE_PC_DEBUG] After table.compute: pc=" << pc << ", dpc_dS=" << dpc_dS << std::endl;
      }
    }
  }
  // For scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
  else
  {
    // Precomputed values should be provided by caller (from computeFlux before local_solver)
    // If not provided, fall back to drainage curve to avoid calling computeTrappedCriticalPhaseVolFraction
    if( precomputedScrt < 0.0 )
    {
      // Precomputed values not provided - return drainage curve value
      real64 const arrayIndex = fields::cappres::ModeIndexType::DRAINAGE;
      if( arrayIndex < static_cast< integer >(capPresKernelWrappers.size()))
      {
        pc = capPresKernelWrappers[arrayIndex].compute( &S, &dpc_dS );
      }
      return pc;
    }

    real64 dpci_dS, dpcd_dS;
    real64 const pci = capPresKernelWrappers[fields::cappres::ModeIndexType::IMBIBITION].compute( &S, &dpci_dS );
    real64 const pcd = capPresKernelWrappers[fields::cappres::ModeIndexType::DRAINAGE].compute( &S, &dpcd_dS );

    real64 const E = killoughCurvatureParam;

    if( isWettingPhase )
    {
      real64 const Smin = wettingCurve.oppositeBoundPhaseVolFraction;
      real64 const Smax = wettingCurve.drainageExtremaPhaseVolFraction;

      if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
      {
        // DRAINAGE_TO_IMBIBITION: transitioning from drainage (low S) to imbibition (high S)
        // Shy should be the minimum historical saturation (where drainage started)
        real64 const Shy = (precomputedShy >= 0.0) ? precomputedShy :
                           ((phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin);

        // Use precomputed value - should have been computed before local_solver
        real64 const Scrt = precomputedScrt;
        // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
        // This matches the forward compute: Swma = 1 - (1 - Scrt) = Scrt
        real64 const Swma = 1 - (1 - Scrt);
        real64 denomF = precomputedDenomF;
        if( LvArray::math::abs( precomputedDenomF ) < 1e-15 )
        {
          denomF = (1. / (Swma - Shy + E) - 1. / E);
        }
        // Guard against division by zero
        real64 F = 0.0;
        if( LvArray::math::abs( denomF ) >= 1e-15 )
        {
          real64 const F_num = (1. / (S - Shy + E) - 1. / E);
          F = F_num / denomF;
        }
        F = LvArray::math::max( F, 0.0 );
        F = LvArray::math::min( F, 1.0 );

        pc = pcd + F * (pci - pcd);
        dpc_dS = dpcd_dS + F * (dpci_dS - dpcd_dS);
      }
      else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
      {
        // IMBIBITION_TO_DRAINAGE: transitioning from imbibition (high S) to drainage (low S)
        // Shy should be the maximum historical saturation (where imbibition started)
        real64 const Shy = (precomputedShy >= 0.0) ? precomputedShy :
                           ((phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax);

        // Use precomputed value - should have been computed before local_solver
        real64 const Scrt = precomputedScrt;

        real64 denomF = precomputedDenomF;
        if( LvArray::math::abs( precomputedDenomF ) < 1e-15 )
        {
          denomF = (1. / (Shy - Scrt + E) - 1. / E);
        }
        // Guard against division by zero
        real64 F = 0.0;
        if( LvArray::math::abs( denomF ) >= 1e-15 )
        {
          real64 const F_num = (1. / (Shy - S + E) - 1. / E);
          F = F_num / denomF;
        }
        F = LvArray::math::max( F, 0.0 );
        F = LvArray::math::min( F, 1.0 );

        pc = pci + F * (pcd - pci);
        dpc_dS = dpci_dS + F * (dpcd_dS - dpci_dS);
      }
      else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
      {
        // IMBIBITION_TO_DRAINAGE_FROM_SCANNING: Secondary drainage scanning curve (Mode 4)
        // Based on Killough (1976) - solve quadratic for ghost departure point
        real64 const H = (phaseMinHistoricalVolFraction > Smin) ? phaseMinHistoricalVolFraction : Smin;
        real64 S_star = phaseMode2PeakVolFraction;
        if( S_star < 1e-12 || S_star >= Smax - 1e-12 )
        {
          if( phaseMaxHistoricalVolFraction > H + 1e-12 && phaseMaxHistoricalVolFraction < Smax - 1e-12 )
          {
            S_star = phaseMaxHistoricalVolFraction;
          }
          else
          {
            S_star = LvArray::math::max( H + 0.1, LvArray::math::min( S, Smax - 0.01 ));
          }
        }

        // Compute F_known from Mode 2 at S_star
        real64 const Scrt = precomputedScrt;
        // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
        real64 const Swma = 1 - (1 - Scrt);
        real64 F_known = 0.0;
        if( S_star <= Swma + 1e-12 )
        {
          real64 const F_num = (1. / (S_star - H + E) - 1. / E);
          real64 const F_denom = (1. / (Swma - H + E) - 1. / E);
          if( LvArray::math::abs( F_denom ) > 1e-12 )
          {
            F_known = F_num / F_denom;
          }
        }
        F_known = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_known ));

        // F_star_target = 1 - F_known
        real64 const F_star_target = 1.0 - F_known;

        // Solve quadratic for x (ghost departure point)
        real64 const a = S_star;
        real64 const b = H;                  // Use H so curve rejoins drainage at H
        real64 const T = F_star_target;
        real64 const coeff_A = (1.0 - T);
        real64 const coeff_B = -(1.0 - T) * (a + b - 2.0 * E) - (1.0 - T) * E;
        real64 const coeff_C = (1.0 - T) * (a - E) * (b - E)
                               - E * E * (1.0 + T)
                               - E * (T * a - b);

        real64 x = S_star;
        real64 discriminant = coeff_B * coeff_B - 4.0 * coeff_A * coeff_C;
        if( discriminant >= 0.0 && LvArray::math::abs( coeff_A ) > 1e-14 )
        {
          real64 sqrt_disc = LvArray::math::sqrt( discriminant );
          real64 x1 = (-coeff_B + sqrt_disc) / (2.0 * coeff_A);
          real64 x2 = (-coeff_B - sqrt_disc) / (2.0 * coeff_A);
          if( x1 > S_star - 1e-8 && x2 > S_star - 1e-8 )
          {
            x = LvArray::math::min( x1, x2 );
          }
          else if( x1 > S_star - 1e-8 )
          {
            x = x1;
          }
          else if( x2 > S_star - 1e-8 )
          {
            x = x2;
          }
        }

        // Compute F_star: F_star = [1/(x - S + E) - 1/E] / [1/(x - H + E) - 1/E]
        real64 F_star = 1.0;
        real64 const F_star_num = (1. / (x - S + E) - 1. / E);
        real64 const F_star_denom = (1. / (x - H + E) - 1. / E);
        if( LvArray::math::abs( F_star_denom ) > 1e-12 )
        {
          F_star = F_star_num / F_star_denom;
          F_star = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_star ));
        }

        // Pc = Pc_Im + F_star * (Pc_Dr - Pc_Im)
        pc = pci + F_star * (pcd - pci);
        dpc_dS = dpci_dS + F_star * (dpcd_dS - dpci_dS);
      }
    }
    else
    {
      // Non-wetting phase
      real64 const Smax = nonWettingCurve.oppositeBoundPhaseVolFraction;
      real64 const Shy = (precomputedShy >= 0.0) ? precomputedShy :
                         ((phaseMaxHistoricalVolFraction < Smax) ? phaseMaxHistoricalVolFraction : Smax);

      if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
      {
        // Use precomputed value - should have been computed before local_solver
        real64 const Scrt = precomputedScrt;

        real64 denomF = precomputedDenomF;
        if( LvArray::math::abs( precomputedDenomF ) < 1e-15 )
        {
          denomF = (1. / (Shy - Scrt + E) - 1. / E);
        }
        // Guard against division by zero
        real64 F = 0.0;
        if( LvArray::math::abs( denomF ) >= 1e-15 )
        {
          real64 const F_num = (1. / (Shy - S + E) - 1. / E);
          F = F_num / denomF;
        }
        F = LvArray::math::max( F, 0.0 );
        F = LvArray::math::min( F, 1.0 );

        pc = pcd + F * (pci - pcd);
        dpc_dS = dpci_dS + F * (dpci_dS - dpcd_dS);
        // Note: there's also a dF/dS term in the non-wetting case, but for inversion we approximate
      }
      else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
      {
        // Use precomputed value - should have been computed before local_solver
        real64 const Scrt = precomputedScrt;
        real64 const Sgma = 1. - Scrt - phaseIntermediateMinVolFraction;

        real64 denomF = precomputedDenomF;
        if( LvArray::math::abs( precomputedDenomF ) < 1e-15 )
        {
          denomF = (1. / (Sgma - Shy + E) - 1. / E);
        }
        // Guard against division by zero
        real64 F = 0.0;
        if( LvArray::math::abs( denomF ) >= 1e-15 )
        {
          real64 const F_num = (1. / (S - Shy + E) - 1. / E);
          F = F_num / denomF;
        }
        F = LvArray::math::max( F, 0.0 );
        F = LvArray::math::min( F, 1.0 );

        pc = pci + F * (pcd - pci);
        dpc_dS = dpcd_dS + F * (dpcd_dS - dpci_dS);
        // Note: there's also a dF/dS term in the non-wetting case, but for inversion we approximate
      }
    }
  }

  return pc;
}

void
TableCapillaryPressureHysteresis::KernelWrapper::computeInv(
  arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFraction,
  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMode2PeakVolFraction,
  arraySlice1d< real64 const, cappres::USD_CAPPRES - 2 > const & phaseTrappedVolFrac,
  arraySlice1d< real64 const, cappres::USD_CAPPRES - 2 > const & phaseCapPressure,
  arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac,
  fields::cappres::ModeIndexType const & mode ) const
{
  GEOS_UNUSED_VAR( phaseTrappedVolFrac );
  LvArray::forValuesInSlice( dPhaseCapPressure_dPhaseVolFrac, []( real64 & val ) { val = 0.0; } );

  constexpr bool ENABLE_COMPUTEINV_PATH_DEBUG = false;
  if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
    std::cout << "[COMPUTEINV_PATH] Entry to computeInv(), mode="
              << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" :
        mode == fields::cappres::ModeIndexType::IMBIBITION ? "IMBIBITION" : "SCANNING") << std::endl;
  }

  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = (PT::WATER < m_phaseOrder.size()) ? m_phaseOrder[PT::WATER] : -1;
  integer const ipOil = (PT::OIL < m_phaseOrder.size()) ? m_phaseOrder[PT::OIL] : -1;
  integer const ipGas = (PT::GAS < m_phaseOrder.size()) ? m_phaseOrder[PT::GAS] : -1;

  if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
    std::cout << "[COMPUTEINV_PATH] Phase indices: ipWater=" << ipWater
              << ", ipOil=" << ipOil << ", ipGas=" << ipGas << std::endl;
    std::cout << "[COMPUTEINV_PATH] After phase indices calculation, continuing..." << std::endl;
  }

  constexpr real64 tol = 1e-9;
  constexpr integer maxIter = 20;
  constexpr real64 minS = 0.0;
  constexpr real64 maxS = 1.0;

  // Determine which phase pairs need inversion
  if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
    std::cout << "[COMPUTEINV_PATH] Checking phase combination..." << std::endl;
    std::cout << "[COMPUTEINV_PATH] ipWater=" << ipWater << ", ipOil=" << ipOil << ", ipGas=" << ipGas << std::endl;
    bool isThreePhase = (ipWater >= 0 && ipOil >= 0 && ipGas >= 0);
    std::cout << "[COMPUTEINV_PATH] Three-phase condition: " << (isThreePhase ? "TRUE" : "FALSE") << std::endl;
    std::cout << "[COMPUTEINV_PATH] About to check three-phase if statement..." << std::endl;
  }
  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
      std::cout << "[COMPUTEINV_PATH] Taking THREE-PHASE path" << std::endl;
    }
    // Three-phase: invert wetting and non-wetting capillary pressures

    // 1. Invert wetting phase (water-oil capillary pressure)
    constexpr real64 pcEpsilon = 1e-10;
    constexpr bool ENABLE_INVERSE_TABLE_DEBUG = false;
    if( ipWater >= 0 && LvArray::math::abs( phaseCapPressure[ipWater] ) > pcEpsilon )
    {
      // For pure DRAINAGE or IMBIBITION modes, use direct table lookup (analytical inverse)
      if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
      {
        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? 0 : 1;
        array1d< real64 > input( 1 );
        input[0] = phaseCapPressure[ipWater];
        auto inputSlice = input.toSliceConst();

        if constexpr (ENABLE_INVERSE_TABLE_DEBUG) {
          real64 S_before = phaseVolFraction[ipWater];
          real64 dS_dPc_before = dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater];
          std::cout << "[INVERSE_TABLE_DEBUG] Cell1 (water), mode=" << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                    << ", arrayIndex=" << arrayIndex << std::endl;
          std::cout << "  Input Pc=" << input[0] << std::endl;
          std::cout << "  S before inverse table=" << S_before << std::endl;
          std::cout << "  dS_dPc before inverse table=" << dS_dPc_before << std::endl;
        }

        // Match TableCapillaryPressure exactly: pass derivative pointer directly to inverse table
        // Note: The inverse table returns dS/dPc, but we store it directly like TableCapillaryPressure does
        // The local solver will call compute() afterwards to get the correct dPc/dS derivative
        phaseVolFraction[ipWater] = m_inverseWettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute(
          inputSlice, &dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] );

        if constexpr (ENABLE_INVERSE_TABLE_DEBUG) {
          real64 S_after = phaseVolFraction[ipWater];
          real64 dS_dPc_after = dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater];
          std::cout << "  S after inverse table=" << S_after << std::endl;
          std::cout << "  dS_dPc after inverse table=" << dS_dPc_after << std::endl;
        }

        phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipWater] - phaseVolFraction[ipGas];
        // Direct lookup complete, skip Newton-Raphson
      }
      else
      {
        // For scanning curves, use analytical inverse if possible, otherwise Newton-Raphson
        // TEMPORARY: Flag to enable analytical inverse for scanning curves
        constexpr bool USE_ANALYTICAL_INVERSE_SCANNING = true;

        real64 S_guess = phaseMaxHistoricalVolFraction[ipWater];
        if( S_guess <= phaseMinHistoricalVolFraction[ipWater] || S_guess< minS || S_guess > maxS )
        {
          // Fall back to drainage/imbibition curve evaluation if available
          real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
          real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
          S_guess = 0.5 * (Smin + Smax);
          S_guess = LvArray::math::max( S_guess, minS );
          S_guess = LvArray::math::min( S_guess, maxS );
        }

        // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
        real64 precomputedScrt_water = -1.0;
        real64 precomputedDenomF_water = 0.0;
        real64 precomputedShy_water = -1.0;

        // Precompute for scanning curves (Mode 2, Mode 3, or Mode 4)
        if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
        {
          real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
          real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
          real64 const E = m_killoughCurvatureParamCapPres;

          if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
          {
            // DRAINAGE_TO_IMBIBITION: Shy should be minimum historical (where drainage started)
            precomputedShy_water = (phaseMinHistoricalVolFraction[ipWater] > Smin) ?
                                   phaseMinHistoricalVolFraction[ipWater] : Smin;

            // Compute Scrt for DRAINAGE_TO_IMBIBITION scanning curve (wetting phase)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_water,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
            real64 const Swma = 1. - (1. - Scrt);
            // Compute denomF for DRAINAGE_TO_IMBIBITION: denomF = (1. / (Swma - Shy + E) - 1. / E)
            precomputedDenomF_water = (1. / (Swma - precomputedShy_water + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
          {
            // IMBIBITION_TO_DRAINAGE: Shy should be maximum historical (where imbibition started)
            precomputedShy_water = (phaseMaxHistoricalVolFraction[ipWater] < Smax) ?
                                   phaseMaxHistoricalVolFraction[ipWater] : Smax;

            // Compute Scrt for IMBIBITION_TO_DRAINAGE scanning curve (wetting phase)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_water,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // Compute denomF for IMBIBITION_TO_DRAINAGE: denomF = (1. / (Shy - Scrt + E) - 1. / E)
            precomputedDenomF_water = (1. / (precomputedShy_water - Scrt + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
          {
            // Mode 4: IMBIBITION_TO_DRAINAGE_FROM_SCANNING
            // H (first reversal point) should be minimum historical (where imbibition started)
            // This is used to compute Scrt for Mode 4
            real64 const H = (phaseMinHistoricalVolFraction[ipWater] > Smin) ?
                             phaseMinHistoricalVolFraction[ipWater] : Smin;

            // Compute Scrt for Mode 4 based on H (first reversal point)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, H,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // For Mode 4, precomputedShy is not used in the same way, but set it to H for consistency
            precomputedShy_water = H;
            // denomF is not used for Mode 4 (it solves quadratic instead), but set to 0 to indicate it's not applicable
            precomputedDenomF_water = 0.0;
          }
        }

        real64 S = S_guess;

        // Try analytical inverse for scanning curves if enabled
        if( USE_ANALYTICAL_INVERSE_SCANNING )
        {
          // Use inverse tables to get bounds: S_dr(Pc) and S_im(Pc)
          array1d< real64 > input_dr( 1 ), input_im( 1 );
          input_dr[0] = phaseCapPressure[ipWater];
          input_im[0] = phaseCapPressure[ipWater];
          auto inputSlice_dr = input_dr.toSliceConst();
          auto inputSlice_im = input_im.toSliceConst();

          real64 dS_dPc_dr = 0.0, dS_dPc_im = 0.0;
          real64 S_dr = m_inverseWettingIntermediateCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE].compute(
            inputSlice_dr, &dS_dPc_dr );
          real64 S_im = m_inverseWettingIntermediateCapillaryPressureKernelWrappers[ModeIndexType::IMBIBITION].compute(
            inputSlice_im, &dS_dPc_im );

          // Clamp to valid range
          S_dr = LvArray::math::max( minS, LvArray::math::min( maxS, S_dr ));
          S_im = LvArray::math::max( minS, LvArray::math::min( maxS, S_im ));

          // For scanning curves, S is between S_dr and S_im (order depends on mode)
          real64 S_low = LvArray::math::min( S_dr, S_im );
          real64 S_high = LvArray::math::max( S_dr, S_im );

          // Use analytical inverse based on F(S) formula
          real64 const E = m_killoughCurvatureParamCapPres;
          real64 const Pc_target = phaseCapPressure[ipWater];

          if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
          {
            // DRAINAGE_TO_IMBIBITION: Pc = Pc_dr(S) + F(S) * [Pc_im(S) - Pc_dr(S)]
            // F(S) = (1. / (S - Shy + E) - 1. / E) / denomF
            // Rearranging: S = Shy - E + 1 / (F(S) * denomF + 1/E)
            real64 const Shy = precomputedShy_water;
            real64 const denomF = precomputedDenomF_water;

            // Initial guess: interpolate between S_dr and S_im
            S = 0.5 * (S_low + S_high);

            // Fixed-point iteration using analytical F(S) formula
            for( integer iter = 0; iter < 10; ++iter )
            {
              // Compute F(S) from current S
              real64 F_S = (1. / (S - Shy + E) - 1. / E) / denomF;
              F_S = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_S ));

              // Compute Pc(S) to check convergence
              real64 pc_computed = this->computeCapillaryPressureForSaturation(
                S, mode, ipWater,
                phaseMinHistoricalVolFraction[ipWater],
                phaseMaxHistoricalVolFraction[ipWater],
                phaseMode2PeakVolFraction[ipWater],
                m_wettingIntermediateCapillaryPressureKernelWrappers,
                m_wettingCurve,
                m_nonWettingCurve,
                m_landParam[ipWater],
                m_phaseIntermediateMinVolFraction,
                m_killoughCurvatureParamCapPres,
                m_jerauldParam_a,
                m_jerauldParam_b,
                true,                     // isWettingPhase
                precomputedScrt_water,
                precomputedDenomF_water,
                precomputedShy_water );

              if( LvArray::math::abs( pc_computed - Pc_target ) < tol )
              {
                break;
              }

              real64 denom_F = F_S * denomF + 1. / E;
              if( LvArray::math::abs( denom_F ) > 1e-12 )
              {
                real64 S_new = Shy - E + 1. / denom_F;
                S_new = LvArray::math::max( S_low, LvArray::math::min( S_high, S_new ));
                S = S_new;
              }
              else
              {
                // Fall back to interpolation if denominator is too small
                S = 0.5 * (S + 0.5 * (S_low + S_high));
              }
            }
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
          {
            // IMBIBITION_TO_DRAINAGE: Pc = Pc_im(S) + F(S) * [Pc_dr(S) - Pc_im(S)]
            // F(S) = (1. / (S - Shy + E) - 1. / E) / denomF
            real64 const Shy = precomputedShy_water;
            real64 const denomF = precomputedDenomF_water;

            S = 0.5 * (S_low + S_high);

            for( integer iter = 0; iter < 10; ++iter )
            {
              real64 F_S = (1. / (S - Shy + E) - 1. / E) / denomF;
              F_S = LvArray::math::max( 0.0, LvArray::math::min( 1.0, F_S ));

              real64 pc_computed = this->computeCapillaryPressureForSaturation(
                S, mode, ipWater,
                phaseMinHistoricalVolFraction[ipWater],
                phaseMaxHistoricalVolFraction[ipWater],
                phaseMode2PeakVolFraction[ipWater],
                m_wettingIntermediateCapillaryPressureKernelWrappers,
                m_wettingCurve,
                m_nonWettingCurve,
                m_landParam[ipWater],
                m_phaseIntermediateMinVolFraction,
                m_killoughCurvatureParamCapPres,
                m_jerauldParam_a,
                m_jerauldParam_b,
                true,                     // isWettingPhase
                precomputedScrt_water,
                precomputedDenomF_water,
                precomputedShy_water );

              if( LvArray::math::abs( pc_computed - Pc_target ) < tol )
              {
                break;
              }

              real64 denom_F = F_S * denomF + 1. / E;
              if( LvArray::math::abs( denom_F ) > 1e-12 )
              {
                real64 S_new = Shy - E + 1. / denom_F;
                S_new = LvArray::math::max( S_low, LvArray::math::min( S_high, S_new ));
                S = S_new;
              }
              else
              {
                S = 0.5 * (S + 0.5 * (S_low + S_high));
              }
            }
          }
          // For IMBIBITION_TO_DRAINAGE_FROM_SCANNING, fall through to Newton-Raphson (more complex F_star formula)
        }

        real64 pc_check = this->computeCapillaryPressureForSaturation(
          S, mode, ipWater,
          phaseMinHistoricalVolFraction[ipWater],
          phaseMaxHistoricalVolFraction[ipWater],
          phaseMode2PeakVolFraction[ipWater],
          m_wettingIntermediateCapillaryPressureKernelWrappers,
          m_wettingCurve,
          m_nonWettingCurve,
          m_landParam[ipWater],
          m_phaseIntermediateMinVolFraction,
          m_killoughCurvatureParamCapPres,
          m_jerauldParam_a,
          m_jerauldParam_b,
          true,               // isWettingPhase
          precomputedScrt_water,
          precomputedDenomF_water,
          precomputedShy_water );

        bool needs_newton = (mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING) ||
                            (LvArray::math::abs( pc_check - phaseCapPressure[ipWater] ) >= tol);

        if( needs_newton )
        {
          real64 residual_prev = 1e30;

          for( integer iter = 0; iter < maxIter; ++iter )
          {
            // Compute Pc(S) using the helper function
            real64 pc_computed = this->computeCapillaryPressureForSaturation(
              S, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingIntermediateCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );

            real64 residual = pc_computed - phaseCapPressure[ipWater];

            if( LvArray::math::abs( residual ) < tol )
            {
              break;
            }

            // Compute derivative dPc/dS using adaptive finite difference
            // Use relative step size to handle different scales
            real64 const dS = LvArray::math::max( 1e-8, 1e-6 * LvArray::math::max( LvArray::math::abs( S ), 1e-6 ));
            real64 S_pert = S + dS;
            if( S_pert > maxS )
            {
              S_pert = LvArray::math::max( S - dS, minS );
            }

            real64 const pc_pert = this->computeCapillaryPressureForSaturation(
              S_pert, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingIntermediateCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );

            real64 const dpc_dS = (pc_pert - pc_computed) / (S_pert - S);

            if( LvArray::math::abs( dpc_dS ) > 1e-12 )
            {
              real64 const newtonStep = -residual / dpc_dS;
              real64 const maxStep = 0.1;               // Limit step size to prevent overshooting
              real64 const step = LvArray::math::max( -maxStep, LvArray::math::min( maxStep, newtonStep ));

              real64 S_new = S + step;

              // Line search: if residual increases, reduce step
              if( iter > 0 && LvArray::math::abs( residual ) > LvArray::math::abs( residual_prev ))
              {
                // Backtrack: use smaller step
                S_new = S + 0.5 * step;
              }

              // Clamp to valid range
              S_new = LvArray::math::max( S_new, minS );
              S_new = LvArray::math::min( S_new, maxS );

              if( LvArray::math::abs( S_new - S ) < 1e-10 &&
                  (S_new <= minS + 1e-10 || S_new >= maxS - 1e-10))
              {
                // Use bisection between current guess and boundary
                S_new = 0.5 * (S + (S_new <= minS + 1e-10 ? minS : maxS));
              }

              residual_prev = residual;
              S = S_new;
            }
            else
            {
              // If derivative is zero or very small, use bisection
              real64 const S_low = LvArray::math::max( S - 0.1, minS );
              real64 const S_high = LvArray::math::min( S + 0.1, maxS );
              S = 0.5 * (S_low + S_high);
              S = LvArray::math::max( S, minS );
              S = LvArray::math::min( S, maxS );
            }
          }

          phaseVolFraction[ipWater] = S;
          phaseVolFraction[ipOil] = 1.0 - S - phaseVolFraction[ipGas];

          // Compute derivative at final S
          real64 dpc_dS_final = 0.0;
          if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
          {
            integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
            m_wettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute( &S, &dpc_dS_final );
          }
          else
          {
            // For scanning curves, compute derivative using finite difference
            real64 const dS_final = 1e-6;
            real64 const S_pert_final = LvArray::math::min( S + dS_final, maxS );
            real64 const pc_final = this->computeCapillaryPressureForSaturation(
              S, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingIntermediateCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );
            real64 const pc_pert_final = this->computeCapillaryPressureForSaturation(
              S_pert_final, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingIntermediateCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );
            dpc_dS_final = (pc_pert_final - pc_final) / dS_final;
          }
          dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] = dpc_dS_final;
        }             // End else block for scanning curves
      }
    }             // Close if (ipWater >= 0) block

    // 2. Invert non-wetting phase (gas-oil capillary pressure)
    if( ipGas >= 0 && LvArray::math::abs( phaseCapPressure[ipGas] ) > pcEpsilon )
    {
      // For pure DRAINAGE or IMBIBITION modes, use direct table lookup (analytical inverse)
      if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
      {
        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? 0 : 1;
        array1d< real64 > input( 1 );
        input[0] = phaseCapPressure[ipGas];
        auto inputSlice = input.toSliceConst();

        if constexpr (ENABLE_INVERSE_TABLE_DEBUG) {
          real64 S_before = phaseVolFraction[ipGas];
          real64 dS_dPc_before = dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas];
          std::cout << "[INVERSE_TABLE_DEBUG] Cell2 (gas), mode=" << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                    << ", arrayIndex=" << arrayIndex << std::endl;
          std::cout << "  Input Pc=" << input[0] << std::endl;
          std::cout << "  S before inverse table=" << S_before << std::endl;
          std::cout << "  dS_dPc before inverse table=" << dS_dPc_before << std::endl;
        }

        // Match TableCapillaryPressure exactly: pass derivative pointer directly to inverse table
        // Note: The inverse table returns dS/dPc, but we store it directly like TableCapillaryPressure does
        // The local solver will call compute() afterwards to get the correct dPc/dS derivative
        phaseVolFraction[ipGas] = m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute(
          inputSlice, &dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas] );

        if constexpr (ENABLE_INVERSE_TABLE_DEBUG) {
          real64 S_after = phaseVolFraction[ipGas];
          real64 dS_dPc_after = dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas];
          std::cout << "  S after inverse table=" << S_after << std::endl;
          std::cout << "  dS_dPc after inverse table=" << dS_dPc_after << std::endl;
        }

        if( ipWater >= 0 && ipOil >= 0 )
        {
          phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipWater] - phaseVolFraction[ipGas];
        }
        // Direct lookup complete, skip Newton-Raphson
      }
      else
      {
        // For scanning curves, use Newton-Raphson
        // Note: For non-wetting phase, the capillary pressure is typically negative
        // and the relationship may be inverted (increasing Pc with decreasing S)
        real64 S_guess = phaseMaxHistoricalVolFraction[ipGas];
        if( S_guess <= phaseMinHistoricalVolFraction[ipGas] || S_guess< minS || S_guess > maxS )
        {
          real64 const Smin = m_nonWettingCurve.imbibitionExtremaPhaseVolFraction;
          real64 const Smax = m_nonWettingCurve.drainageExtremaPhaseVolFraction;
          S_guess = 0.5 * (Smin + Smax);
          S_guess = LvArray::math::max( S_guess, minS );
          S_guess = LvArray::math::min( S_guess, maxS );
        }

        // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
        real64 precomputedScrt_gas = -1.0;
        real64 precomputedDenomF_gas = 0.0;
        real64 precomputedShy_gas = -1.0;

        // Only precompute for scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
        if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
        {
          real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
          precomputedShy_gas = (phaseMaxHistoricalVolFraction[ipGas] < Smax) ?
                               phaseMaxHistoricalVolFraction[ipGas] : Smax;
          real64 const E = m_killoughCurvatureParamCapPres;

          if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
          {
            // Compute Scrt for DRAINAGE_TO_IMBIBITION scanning curve
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, precomputedShy_gas,
                                                                        m_landParam[ipGas],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_gas = Scrt;
            // Compute denomF for DRAINAGE_TO_IMBIBITION: denomF = (1. / (Shy - Scrt + E) - 1. / E)
            precomputedDenomF_gas = (1. / (precomputedShy_gas - Scrt + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
          {
            // Compute Scrt for IMBIBITION_TO_DRAINAGE scanning curve
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_gas,
                                                                        m_landParam[ipGas],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_gas = Scrt;
            real64 const Sgma = 1. - Scrt - m_phaseIntermediateMinVolFraction;
            // Compute denomF for IMBIBITION_TO_DRAINAGE: denomF = (1. / (Sgma - Shy + E) - 1. / E)
            precomputedDenomF_gas = (1. / (Sgma - precomputedShy_gas + E) - 1. / E);
          }
        }

        real64 S = S_guess;
        real64 residual_prev = 1e30;

        for( integer iter = 0; iter < maxIter; ++iter )
        {
          real64 pc_computed = this->computeCapillaryPressureForSaturation(
            S, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );

          // For non-wetting phase, capillary pressure is typically multiplied by -1
          // Check which sign to use based on the input
          real64 residual = pc_computed - phaseCapPressure[ipGas];

          if( LvArray::math::abs( residual ) < tol )
          {
            break;
          }

          // Compute derivative dPc/dS using adaptive finite difference
          real64 const dS = LvArray::math::max( 1e-8, 1e-6 * LvArray::math::max( LvArray::math::abs( S ), 1e-6 ));
          real64 S_pert = S + dS;
          if( S_pert > maxS )
          {
            S_pert = LvArray::math::max( S - dS, minS );
          }

          real64 const pc_pert = this->computeCapillaryPressureForSaturation(
            S_pert, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );

          real64 const dpc_dS = (pc_pert - pc_computed) / (S_pert - S);

          if( LvArray::math::abs( dpc_dS ) > 1e-12 )
          {
            real64 const newtonStep = -residual / dpc_dS;
            real64 const maxStep = 0.1;                 // Limit step size to prevent overshooting
            real64 const step = LvArray::math::max( -maxStep, LvArray::math::min( maxStep, newtonStep ));

            real64 S_new = S + step;

            // Line search: if residual increases, reduce step
            if( iter > 0 && LvArray::math::abs( residual ) > LvArray::math::abs( residual_prev ))
            {
              // Backtrack: use smaller step
              S_new = S + 0.5 * step;
            }

            // Clamp to valid range
            S_new = LvArray::math::max( S_new, minS );
            S_new = LvArray::math::min( S_new, maxS );

            if( LvArray::math::abs( S_new - S ) < 1e-10 &&
                (S_new <= minS + 1e-10 || S_new >= maxS - 1e-10))
            {
              // Use bisection between current guess and boundary
              S_new = 0.5 * (S + (S_new <= minS + 1e-10 ? minS : maxS));
            }

            residual_prev = residual;
            S = S_new;
          }
          else
          {
            // If derivative is zero or very small, use bisection
            real64 const S_low = LvArray::math::max( S - 0.1, minS );
            real64 const S_high = LvArray::math::min( S + 0.1, maxS );
            S = 0.5 * (S_low + S_high);
            S = LvArray::math::max( S, minS );
            S = LvArray::math::min( S, maxS );
          }
        }

        phaseVolFraction[ipGas] = S;

        if( ipWater >= 0 && ipOil >= 0 )
        {
          phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipWater] - phaseVolFraction[ipGas];
          phaseVolFraction[ipOil] = LvArray::math::max( phaseVolFraction[ipOil], 0.0 );
          phaseVolFraction[ipOil] = LvArray::math::min( phaseVolFraction[ipOil], 1.0 );
        }

        // Compute derivative at final S
        real64 dpc_dS_final = 0.0;
        if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
        {
          integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? fields::cappres::ModeIndexType::DRAINAGE : fields::cappres::ModeIndexType::IMBIBITION;
          m_nonWettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute( &S, &dpc_dS_final );
        }
        else
        {
          // For scanning curves, compute derivative using finite difference
          real64 const dS_final = 1e-6;
          real64 const S_pert_final = LvArray::math::min( S + dS_final, maxS );
          real64 const pc_final = this->computeCapillaryPressureForSaturation(
            S, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );
          real64 const pc_pert_final = this->computeCapillaryPressureForSaturation(
            S_pert_final, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_nonWettingIntermediateCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );
          dpc_dS_final = (pc_pert_final - pc_final) / dS_final;
        }
        dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas] = dpc_dS_final;
      }               // End else block for scanning curves
    }
  }
  // CRITICAL: Add debug right after three-phase block to see if we reach here
  if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
    std::cout << "[COMPUTEINV_PATH] *** AFTER three-phase block, before else-if checks ***" << std::endl;
    std::cout << "[COMPUTEINV_PATH] ipWater=" << ipWater << ", ipOil=" << ipOil << ", ipGas=" << ipGas << std::endl;
    std::cout << "[COMPUTEINV_PATH] Three-phase block was SKIPPED (condition was false)" << std::endl;
    std::cout << "[COMPUTEINV_PATH] About to check if (ipWater < 0)..." << std::endl;
  }
  if( ipWater < 0 )
  {
    if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
      std::cout << "[COMPUTEINV_PATH] *** TAKING ipWater < 0 path (oil-gas) ***" << std::endl;
      std::cout << "[COMPUTEINV_PATH] This should NOT happen with ipWater=0!" << std::endl;
    }
    if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
      std::cout << "[COMPUTEINV_PATH] After three-phase block, checking else-if conditions..." << std::endl;
      std::cout << "[COMPUTEINV_PATH] ipWater=" << ipWater << ", condition (ipWater < 0)=" << (ipWater < 0 ? "TRUE" : "FALSE") << std::endl;
      std::cout << "[COMPUTEINV_PATH] *** TAKING ipWater < 0 path (oil-gas) ***" << std::endl;
      std::cout << "[COMPUTEINV_PATH] This should NOT happen with ipWater=0!" << std::endl;
    }
    if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
      std::cout << "[COMPUTEINV_PATH] *** TAKING ipWater < 0 path (oil-gas) ***" << std::endl;
      std::cout << "[COMPUTEINV_PATH] This should NOT happen with ipWater=0!" << std::endl;
    }
    // Two-phase: oil-gas (non-wetting phase)
    // Similar to above but simpler
    constexpr real64 pcEpsilon_gas = 1e-10;
    if( ipGas >= 0 && LvArray::math::abs( phaseCapPressure[ipGas] ) > pcEpsilon_gas )
    {
      // For pure DRAINAGE or IMBIBITION modes, use direct table lookup (analytical inverse)
      if( mode == fields::cappres::ModeIndexType::DRAINAGE || mode == fields::cappres::ModeIndexType::IMBIBITION )
      {
        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? 0 : 1;
        array1d< real64 > input( 1 );
        input[0] = phaseCapPressure[ipGas];
        auto inputSlice = input.toSliceConst();
        // Match TableCapillaryPressure: no clamping, let the table handle bounds
        phaseVolFraction[ipGas] = m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers[arrayIndex].compute(
          inputSlice, &dPhaseCapPressure_dPhaseVolFrac[ipGas][ipGas] );
        if( ipOil >= 0 )
        {
          phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipGas];
        }
        // Direct lookup complete, skip Newton-Raphson
      }
      else
      {
        // For scanning curves, use Newton-Raphson
        real64 S_guess = phaseMaxHistoricalVolFraction[ipGas];
        if( S_guess < minS || S_guess > maxS )
        {
          S_guess = 0.5;
        }

        // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
        real64 precomputedScrt_gas = -1.0;
        real64 precomputedDenomF_gas = 0.0;
        real64 precomputedShy_gas = -1.0;

        // Only precompute for scanning curves (DRAINAGE_TO_IMBIBITION or IMBIBITION_TO_DRAINAGE)
        if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
        {
          real64 const Smax = m_nonWettingCurve.oppositeBoundPhaseVolFraction;
          precomputedShy_gas = (phaseMaxHistoricalVolFraction[ipGas] < Smax) ?
                               phaseMaxHistoricalVolFraction[ipGas] : Smax;
          real64 const E = m_killoughCurvatureParamCapPres;

          if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
          {
            // Compute Scrt for DRAINAGE_TO_IMBIBITION scanning curve
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_nonWettingCurve, precomputedShy_gas,
                                                                        m_landParam[ipGas],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_gas = Scrt;
            // Compute denomF for DRAINAGE_TO_IMBIBITION: denomF = (1. / (Shy - Scrt + E) - 1. / E)
            precomputedDenomF_gas = (1. / (precomputedShy_gas - Scrt + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
          {
            // Compute Scrt for IMBIBITION_TO_DRAINAGE scanning curve
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_gas,
                                                                        m_landParam[ipGas],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_gas = Scrt;
            real64 const Sgma = 1. - Scrt - m_phaseIntermediateMinVolFraction;
            // Compute denomF for IMBIBITION_TO_DRAINAGE: denomF = (1. / (Sgma - Shy + E) - 1. / E)
            precomputedDenomF_gas = (1. / (Sgma - precomputedShy_gas + E) - 1. / E);
          }
        }

        real64 S = S_guess;
        for( integer iter = 0; iter < maxIter; ++iter )
        {
          real64 pc_computed = this->computeCapillaryPressureForSaturation(
            S, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_wettingNonWettingCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );

          real64 residual = pc_computed - phaseCapPressure[ipGas];

          if( LvArray::math::abs( residual ) < tol )
          {
            break;
          }

          real64 const dS = 1e-6;
          real64 const S_pert = LvArray::math::min( S + dS, maxS );
          real64 const pc_pert = this->computeCapillaryPressureForSaturation(
            S_pert, mode, ipGas,
            phaseMinHistoricalVolFraction[ipGas],
            phaseMaxHistoricalVolFraction[ipGas],
            phaseMode2PeakVolFraction[ipGas],
            m_wettingNonWettingCapillaryPressureKernelWrappers,
            m_wettingCurve,
            m_nonWettingCurve,
            m_landParam[ipGas],
            m_phaseIntermediateMinVolFraction,
            m_killoughCurvatureParamCapPres,
            m_jerauldParam_a,
            m_jerauldParam_b,
            false,                 // isWettingPhase = false
            precomputedScrt_gas,
            precomputedDenomF_gas,
            precomputedShy_gas );

          real64 const dpc_dS = (pc_pert - pc_computed) / dS;

          if( LvArray::math::abs( dpc_dS ) > 1e-12 )
          {
            S = S - residual / dpc_dS;
            S = LvArray::math::max( S, minS );
            S = LvArray::math::min( S, maxS );
          }
          else
          {
            S = S - residual * 1e-6;
            S = LvArray::math::max( S, minS );
            S = LvArray::math::min( S, maxS );
          }
        }

        phaseVolFraction[ipGas] = S;
        if( ipOil >= 0 )
        {
          phaseVolFraction[ipOil] = 1.0 - S;
        }
      }               // End else block for scanning curves
    }
  }
  else
  {
    if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
      std::cout << "[COMPUTEINV_PATH] *** REACHED TWO-PHASE path (else block) ***" << std::endl;
      std::cout << "[COMPUTEINV_PATH] ipWater=" << ipWater << ", ipOil=" << ipOil << ", ipGas=" << ipGas << std::endl;
      std::cout << "[COMPUTEINV_PATH] About to check if (ipWater >= 0)" << std::endl;
      std::cout << "[COMPUTEINV_PATH] This is the path that should work for Mode 0 (DRAINAGE)!" << std::endl;
    }
    // Two-phase: water-oil or water-gas (wetting phase)
    // Copy exact pattern from TableCapillaryPressure::computeInv
    if( ipWater >= 0 )
    {
      if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
        std::cout << "[COMPUTEINV_PATH] Two-phase: ipWater >= 0, entering water phase block" << std::endl;
        std::cout << "[COMPUTEINV_PATH] About to check mode: mode=" << static_cast< int >(mode)
                  << ", DRAINAGE=" << static_cast< int >(fields::cappres::ModeIndexType::DRAINAGE)
                  << ", IMBIBITION=" << static_cast< int >(fields::cappres::ModeIndexType::IMBIBITION) << std::endl;
      }
      // For pure DRAINAGE or IMBIBITION modes, use direct table lookup (analytical inverse)
      // CRITICAL DEBUG: Check if mode condition is actually true
      bool isDrainage = (mode == fields::cappres::ModeIndexType::DRAINAGE);
      bool isImbibition = (mode == fields::cappres::ModeIndexType::IMBIBITION);
      bool isPureMode = (isDrainage || isImbibition);
      if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
        std::cout << "[COMPUTEINV_PATH] Mode condition check:" << std::endl;
        std::cout << "[COMPUTEINV_PATH]   mode (as int)=" << static_cast< int >(mode) << std::endl;
        std::cout << "[COMPUTEINV_PATH]   DRAINAGE (as int)=" << static_cast< int >(fields::cappres::ModeIndexType::DRAINAGE) << std::endl;
        std::cout << "[COMPUTEINV_PATH]   IMBIBITION (as int)=" << static_cast< int >(fields::cappres::ModeIndexType::IMBIBITION) << std::endl;
        std::cout << "[COMPUTEINV_PATH]   isDrainage=" << (isDrainage ? "TRUE" : "FALSE") << std::endl;
        std::cout << "[COMPUTEINV_PATH]   isImbibition=" << (isImbibition ? "TRUE" : "FALSE") << std::endl;
        std::cout << "[COMPUTEINV_PATH]   isPureMode (DRAINAGE || IMBIBITION)=" << (isPureMode ? "TRUE" : "FALSE") << std::endl;
      }
      if( isPureMode )
      {
        if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
          std::cout << "[COMPUTEINV_PATH] Using pure DRAINAGE/IMBIBITION path (same as working case)" << std::endl;
          std::cout << "[COMPUTEINV_PATH] mode=" << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION") << std::endl;
        }
        integer const arrayIndex = (mode == fields::cappres::ModeIndexType::DRAINAGE) ? 0 : 1;
        if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
          std::cout << "[COMPUTEINV_PATH] arrayIndex=" << arrayIndex << std::endl;
        }

        // Match TableCapillaryPressure exactly (but skip the forward table call which passes Pc instead of S)
        real64 capPresWater = phaseCapPressure[ipWater];
        if constexpr (ENABLE_COMPUTEINV_PATH_DEBUG) {
          std::cout << "[COMPUTEINV_PATH] Input Pc=" << capPresWater << std::endl;
          std::cout << "[COMPUTEINV_PATH] phaseVolFraction[ipWater] before=" << phaseVolFraction[ipWater] << std::endl;
        }
        array1d< real64 > input( 1 );
        input[0] = capPresWater;
        auto inputSlice = input.toSliceConst();

        constexpr bool ENABLE_COMPUTEINV_DEBUG = false;
        if( ENABLE_COMPUTEINV_DEBUG )
        {
          std::cout << "[COMPUTEINV_DEBUG] Two-phase, mode=" << (mode == fields::cappres::ModeIndexType::DRAINAGE ? "DRAINAGE" : "IMBIBITION")
                    << ", arrayIndex=" << arrayIndex << std::endl;
          std::cout << "[COMPUTEINV_DEBUG] Input Pc=" << capPresWater << std::endl;
          std::cout << "[COMPUTEINV_DEBUG] Derivative before call=" << dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] << std::endl;
        }

        // Use only the inverse table call (the forward table call in TableCapillaryPressure passes Pc instead of S, which is incorrect)
        // The derivative array is already zeroed at the start of computeInv
        // real64 dS_dPc_before = dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater];
        real64 S_before_inv = phaseVolFraction[ipWater];
        if( ENABLE_COMPUTEINV_DEBUG )
        {
          std::cout << "[COMPUTEINV_DEBUG] Two-phase DRAINAGE/IMBIBITION: Before inverse table call" << std::endl;
          std::cout << "  S_before=" << S_before_inv << std::endl;
          std::cout << "  Pc_input=" << capPresWater << std::endl;
          std::cout << "  arrayIndex=" << arrayIndex << std::endl;
          std::cout << "  inputSlice[0]=" << inputSlice[0] << std::endl;
          std::cout << "  wrapper.size()=" << m_inverseWettingNonWettingCapillaryPressureKernelWrappers.size() << std::endl;
        }
        real64 S_returned = m_inverseWettingNonWettingCapillaryPressureKernelWrappers[arrayIndex].compute(
          inputSlice,
          &(dPhaseCapPressure_dPhaseVolFrac)[ipWater][ipWater] );
        phaseVolFraction[ipWater] = S_returned;
        real64 dS_dPc_after = dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater];
        real64 S_after_inv = phaseVolFraction[ipWater];

        if( ENABLE_COMPUTEINV_DEBUG )
        {
          std::cout << "[COMPUTEINV_DEBUG] After inverse table call:" << std::endl;
          std::cout << "  S_returned=" << S_returned << std::endl;
          std::cout << "  S_after=" << S_after_inv << std::endl;
          std::cout << "  dS_dPc_after=" << dS_dPc_after << std::endl;
          if( std::abs( S_after_inv ) < 1e-10 && capPresWater > 3000 && capPresWater < 100000 )
          {
            std::cout << "[COMPUTEINV_DEBUG] ERROR: Pc=" << capPresWater
                      << " is within valid range but inverse table returned S=0!" << std::endl;
            std::cout << "[COMPUTEINV_DEBUG] This will cause zero Jacobian and Pc_int won't update!" << std::endl;
          }
        }

        // Set non-wetting phase - match TableCapillaryPressure
        if( ipGas >= 0 )
        {
          phaseVolFraction[ipGas] = 1.0 - phaseVolFraction[ipWater];
        }
        else if( ipOil >= 0 )
        {
          phaseVolFraction[ipOil] = 1.0 - phaseVolFraction[ipWater];
        }
        // Direct lookup complete, skip Newton-Raphson
      }
      else
      {
        constexpr bool ENABLE_TWOPHASE_SCANNINGINV_DEBUG = false;

        // For scanning curves, use Newton-Raphson
        real64 S_guess = phaseMaxHistoricalVolFraction[ipWater];
        if( S_guess < minS || S_guess > maxS )
        {
          real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
          real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
          S_guess = 0.5 * (Smin + Smax);
          S_guess = LvArray::math::max( S_guess, minS );
          S_guess = LvArray::math::min( S_guess, maxS );
        }
        if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
          std::cout << "[2P_SCAN_INV] Mode=" << static_cast< int >(mode)
                    << ", Pc_target=" << phaseCapPressure[ipWater]
                    << ", S_guess=" << S_guess
                    << ", S_minHist=" << phaseMinHistoricalVolFraction[ipWater]
                    << ", S_maxHist=" << phaseMaxHistoricalVolFraction[ipWater]
                    << ", Smin_curve=" << m_wettingCurve.oppositeBoundPhaseVolFraction
                    << ", Smax_curve=" << m_wettingCurve.drainageExtremaPhaseVolFraction
                    << std::endl;
        }

        // Precompute fixed parameters for scanning curves (these don't change during Newton-Raphson)
        real64 precomputedScrt_water = -1.0;
        real64 precomputedDenomF_water = 0.0;
        real64 precomputedShy_water = -1.0;

        // Precompute for scanning curves (Mode 2, Mode 3, or Mode 4)
        if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE ||
            mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
        {
          real64 const Smin = m_wettingCurve.oppositeBoundPhaseVolFraction;
          real64 const Smax = m_wettingCurve.drainageExtremaPhaseVolFraction;
          real64 const E = m_killoughCurvatureParamCapPres;

          if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
          {
            // DRAINAGE_TO_IMBIBITION: Shy should be minimum historical (where drainage started)
            precomputedShy_water = (phaseMinHistoricalVolFraction[ipWater] > Smin) ?
                                   phaseMinHistoricalVolFraction[ipWater] : Smin;

            // Compute Scrt for DRAINAGE_TO_IMBIBITION scanning curve (wetting phase)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_water,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
            real64 const Swma = 1. - (1. - Scrt);
            // Compute denomF for DRAINAGE_TO_IMBIBITION: denomF = (1. / (Swma - Shy + E) - 1. / E)
            precomputedDenomF_water = (1. / (Swma - precomputedShy_water + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
          {
            // IMBIBITION_TO_DRAINAGE: Shy should be maximum historical (where imbibition started)
            precomputedShy_water = (phaseMaxHistoricalVolFraction[ipWater] < Smax) ?
                                   phaseMaxHistoricalVolFraction[ipWater] : Smax;

            // Compute Scrt for IMBIBITION_TO_DRAINAGE scanning curve (wetting phase)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, precomputedShy_water,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // Compute denomF for IMBIBITION_TO_DRAINAGE: denomF = (1. / (Shy - Scrt + E) - 1. / E)
            precomputedDenomF_water = (1. / (precomputedShy_water - Scrt + E) - 1. / E);
          }
          else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
          {
            // Mode 4: IMBIBITION_TO_DRAINAGE_FROM_SCANNING
            // H (first reversal point) should be minimum historical (where imbibition started)
            // This is used to compute Scrt for Mode 4
            real64 const H = (phaseMinHistoricalVolFraction[ipWater] > Smin) ?
                             phaseMinHistoricalVolFraction[ipWater] : Smin;

            // Compute Scrt for Mode 4 based on H (first reversal point)
            real64 Scrt = 0.0;
            KilloughHysteresis::computeTrappedCriticalPhaseVolFraction( m_wettingCurve, H,
                                                                        m_landParam[ipWater],
                                                                        m_jerauldParam_a, m_jerauldParam_b,
                                                                        Scrt );
            precomputedScrt_water = Scrt;
            // For Mode 4, precomputedShy is not used in the same way, but set it to H for consistency
            precomputedShy_water = H;
            // denomF is not used for Mode 4 (it solves quadratic instead), but set to 0 to indicate it's not applicable
            precomputedDenomF_water = 0.0;
          }
        }

        if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
          std::cout << "[2P_SCAN_INV] Precomputed: Shy=" << precomputedShy_water
                    << ", Scrt=" << precomputedScrt_water
                    << ", denomF=" << precomputedDenomF_water
                    << ", E=" << m_killoughCurvatureParamCapPres
                    << ", Somin=" << m_phaseIntermediateMinVolFraction
                    << std::endl;
          // Evaluate Pc at a few points to see the curve shape
          real64 const S_lo = precomputedShy_water;
          // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
          real64 const Swma = 1. - (1. - precomputedScrt_water);
          std::cout << "[2P_SCAN_INV] Shy=" << S_lo << ", Swma=" << Swma << std::endl;
          for( int probe = 0; probe <= 5; ++probe )
          {
            real64 S_probe = S_lo + (Swma - S_lo) * probe / 5.0;
            S_probe = LvArray::math::max( minS, LvArray::math::min( maxS, S_probe ));
            real64 pc_probe = this->computeCapillaryPressureForSaturation(
              S_probe, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingNonWettingCapillaryPressureKernelWrappers,
              m_wettingCurve, m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a, m_jerauldParam_b,
              true, precomputedScrt_water, precomputedDenomF_water, precomputedShy_water );
            std::cout << "[2P_SCAN_INV]   S=" << S_probe << " -> Pc=" << pc_probe << std::endl;
          }
        }

        real64 S_lo_bracket, S_hi_bracket;
        real64 Swma_debug = -1.0;
        bool degenerateBracket = false;
        if( mode == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION )
        {
          // Mode 2: scanning curve valid in [Shy, Swma]
          // For wetting curve, Scrt IS the max wetting saturation (Swma = Scrt)
          real64 const Swma = 1. - (1. - precomputedScrt_water);
          Swma_debug = Swma;
          S_lo_bracket = precomputedShy_water;
          S_hi_bracket = Swma;
          // Guard: if Swma <= Shy, the scanning curve has no valid range (degenerate case).
          // Fall back to returning S = Shy directly.
          if( Swma <= precomputedShy_water + 1e-12 )
          {
            degenerateBracket = true;
          }
        }
        else if( mode == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE )
        {
          // Mode 3: scanning curve valid in [Scrt_end, Shy]
          S_lo_bracket = precomputedScrt_water;
          S_hi_bracket = precomputedShy_water;
          // Guard: if Shy <= Scrt, degenerate case
          if( precomputedShy_water <= precomputedScrt_water + 1e-12 )
          {
            degenerateBracket = true;
          }
        }
        else
        {
          array1d< real64 > input_inv( 1 );
          input_inv[0] = phaseCapPressure[ipWater];
          auto inputSlice_inv = input_inv.toSliceConst();

          real64 dS_dPc_dr_unused = 0.0, dS_dPc_im_unused = 0.0;
          real64 S_dr = m_inverseWettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE].compute(
            inputSlice_inv, &dS_dPc_dr_unused );
          real64 S_im = m_inverseWettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::IMBIBITION].compute(
            inputSlice_inv, &dS_dPc_im_unused );
          S_dr = LvArray::math::max( minS, LvArray::math::min( maxS, S_dr ));
          S_im = LvArray::math::max( minS, LvArray::math::min( maxS, S_im ));
          S_lo_bracket = LvArray::math::min( S_dr, S_im );
          S_hi_bracket = LvArray::math::max( S_dr, S_im );
        }

        S_lo_bracket = LvArray::math::max( minS, S_lo_bracket );
        S_hi_bracket = LvArray::math::min( maxS, S_hi_bracket );

        // Ensure S_lo <= S_hi (safety: swap if somehow still inverted after clamping)
        if( S_lo_bracket > S_hi_bracket )
        {
          real64 tmp = S_lo_bracket;
          S_lo_bracket = S_hi_bracket;
          S_hi_bracket = tmp;
        }

        if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
          std::cout << "[2P_SCAN_INV] Bracket: S_lo=" << S_lo_bracket
                    << ", S_hi=" << S_hi_bracket
                    << ", Shy=" << precomputedShy_water
                    << ", Swma=" << Swma_debug
                    << ", Scrt=" << precomputedScrt_water
                    << (degenerateBracket ? " DEGENERATE" : "")
                    << std::endl;
        }

        // Handle degenerate bracket: if Swma <= Shy (Mode 2) or Shy <= Scrt (Mode 3),
        // the scanning curve has no valid range. Return S = Shy (the departure point).
        // The Killough formula at S=Shy gives F=0, so Pc = Pc_dr(Shy), which is consistent
        // with the drainage curve at that point.
        real64 const Pc_target = phaseCapPressure[ipWater];

        if( degenerateBracket )
        {
          // Return S = Shy (departure point).
          // The Killough formula at Shy gives F=0, so Pc = Pc_dr(Shy), consistent with drainage.
          phaseVolFraction[ipWater] = precomputedShy_water;
          if( ipGas >= 0 )
          {
            phaseVolFraction[ipGas] = 1.0 - precomputedShy_water;
          }
          else if( ipOil >= 0 )
          {
            phaseVolFraction[ipOil] = 1.0 - precomputedShy_water;
          }
          // Use drainage table derivative at Shy → dS/dPc = 1/(dPc/dS)
          real64 dpc_dS_shy = 0.0;
          m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE].compute(
            &precomputedShy_water, &dpc_dS_shy );
          if( LvArray::math::abs( dpc_dS_shy ) > 1e-14 )
          {
            dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] = 1.0 / dpc_dS_shy;
          }
          else
          {
            dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] = 0.0;
          }
          constexpr bool ENABLE_SCANINV_SUMMARY_DEGEN = false;
          if constexpr (ENABLE_SCANINV_SUMMARY_DEGEN) {
            std::cout << "[SCAN_INV] mode=" << static_cast< int >(mode)
                      << " Pc_tgt=" << Pc_target
                      << " Shy=" << precomputedShy_water
                      << " S=" << precomputedShy_water
                      << " dS/dPc=" << dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater]
                      << " DEGENERATE_BRACKET"
                      << std::endl;
          }
        }

        if( !degenerateBracket )
        {
          // Step 2: Evaluate Pc at bracket endpoints to determine sign of residual

          // Helper macro-like inline evaluation (avoids lambda for GEOS_HOST_DEVICE compatibility)
          #define EVAL_SCANNING_PC( S_eval ) \
            this->computeCapillaryPressureForSaturation( \
              (S_eval), mode, ipWater, \
              phaseMinHistoricalVolFraction[ipWater], \
              phaseMaxHistoricalVolFraction[ipWater], \
              phaseMode2PeakVolFraction[ipWater], \
              m_wettingNonWettingCapillaryPressureKernelWrappers, \
              m_wettingCurve, m_nonWettingCurve, \
              m_landParam[ipWater], \
              m_phaseIntermediateMinVolFraction, \
              m_killoughCurvatureParamCapPres, \
              m_jerauldParam_a, m_jerauldParam_b, \
              true, \
              precomputedScrt_water, \
              precomputedDenomF_water, \
              precomputedShy_water )

          real64 Pc_lo = EVAL_SCANNING_PC( S_lo_bracket );
          real64 Pc_hi = EVAL_SCANNING_PC( S_hi_bracket );
          real64 res_lo = Pc_lo - Pc_target;
          real64 res_hi = Pc_hi - Pc_target;

          if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
            std::cout << "[2P_SCAN_INV] Pc_target=" << Pc_target
                      << ", Pc_lo=" << Pc_lo << " (res=" << res_lo << ")"
                      << ", Pc_hi=" << Pc_hi << " (res=" << res_hi << ")" << std::endl;
          }

          // Early return: if the root is already at a bracket endpoint
          // (i.e., Pc_target ≈ Pc_scan(Shy) or Pc_target ≈ Pc_scan(Swma)),
          // return that endpoint directly. This avoids the bisection sign-test
          // bug where res_lo ≈ 0 causes the update rule `residual * res_lo < 0`
          // to take the wrong branch (0 * negative = 0, NOT < 0), pushing the
          // bracket away from the root and converging to garbage.
          bool skipBisection = false;
          real64 S = 0.5 * (S_lo_bracket + S_hi_bracket);               // default midpoint

          if( LvArray::math::abs( res_lo ) < tol )
          {
            S = S_lo_bracket;
            skipBisection = true;
            if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
              std::cout << "[2P_SCAN_INV] Root at lower bracket endpoint: S=" << S
                        << ", res_lo=" << res_lo << std::endl;
            }
          }
          else if( LvArray::math::abs( res_hi ) < tol )
          {
            S = S_hi_bracket;
            skipBisection = true;
            if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
              std::cout << "[2P_SCAN_INV] Root at upper bracket endpoint: S=" << S
                        << ", res_hi=" << res_hi << std::endl;
            }
          }

          if( !skipBisection )
          {
            // Ensure bracket contains the root (res_lo and res_hi should have opposite signs)
            // If not, Pc_target is at the boundary of the scanning curve.
            // Clamp to the nearest endpoint (don't expand outside valid range).
            if( res_lo * res_hi > 0 )
            {
              // Root not bracketed — Pc_target is at or beyond scanning curve boundary
              // Return the endpoint with the smallest residual
              if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
                std::cout << "[2P_SCAN_INV] Root not bracketed: res_lo=" << res_lo
                          << ", res_hi=" << res_hi << " — using nearest endpoint" << std::endl;
              }
              // Pick the S with smallest |residual| and skip bisection
              if( LvArray::math::abs( res_lo ) < LvArray::math::abs( res_hi ))
              {
                S_hi_bracket = S_lo_bracket;                 // collapse bracket to lo
              }
              else
              {
                S_lo_bracket = S_hi_bracket;                 // collapse bracket to hi
              }
            }

            // Step 3: Bisection (guaranteed convergence within bracket)
            S = 0.5 * (S_lo_bracket + S_hi_bracket);             // start at midpoint
            constexpr integer maxBisectIter = 50;             // 50 bisections give ~15 digits of precision

            for( integer iter = 0; iter < maxBisectIter; ++iter )
            {
              real64 const pc_computed = EVAL_SCANNING_PC( S );
              real64 const residual = pc_computed - Pc_target;

              if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
                if( iter <= 3 || iter == maxBisectIter - 1 || LvArray::math::abs( residual ) < tol )
                {
                  std::cout << "[2P_SCAN_INV] Bisect iter=" << iter
                            << ", S=" << S << ", Pc(S)=" << pc_computed
                            << ", residual=" << residual
                            << ", bracket=[" << S_lo_bracket << "," << S_hi_bracket << "]" << std::endl;
                }
              }

              if( LvArray::math::abs( residual ) < tol )
              {
                if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
                  std::cout << "[2P_SCAN_INV] Converged at iter=" << iter << ", S=" << S << std::endl;
                }
                break;
              }

              // Update bracket — use <= to correctly handle res_lo ≈ 0
              // (when the root is very near S_lo, residual * res_lo ≈ 0
              //  and the < test would take the wrong branch)
              if( residual * res_lo <= 0 )
              {
                S_hi_bracket = S;
                res_hi = residual;
              }
              else
              {
                S_lo_bracket = S;
                res_lo = residual;
              }

              // Bisection step
              S = 0.5 * (S_lo_bracket + S_hi_bracket);

              // Check if bracket is too small
              if( S_hi_bracket - S_lo_bracket < 1e-12 )
              {
                break;
              }
            }
          }               // end if (!skipBisection)

                        #undef EVAL_SCANNING_PC

          if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
            real64 pc_final_check = this->computeCapillaryPressureForSaturation(
              S, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingNonWettingCapillaryPressureKernelWrappers,
              m_wettingCurve, m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a, m_jerauldParam_b,
              true, precomputedScrt_water, precomputedDenomF_water, precomputedShy_water );
            std::cout << "[2P_SCAN_INV] FINAL: S=" << S
                      << ", Pc(S)=" << pc_final_check
                      << ", Pc_target=" << phaseCapPressure[ipWater]
                      << ", error=" << (pc_final_check - phaseCapPressure[ipWater]) << std::endl;
          }

          phaseVolFraction[ipWater] = S;
          // Set non-wetting phase
          if( ipGas >= 0 )
          {
            phaseVolFraction[ipGas] = 1.0 - S;
          }
          else if( ipOil >= 0 )
          {
            phaseVolFraction[ipOil] = 1.0 - S;
          }

          // Compute derivative dPc/dS at final S using finite difference
          // (matching the three-phase scanning curve path)
          {
            real64 const dS_final_requested = 1e-6;
            real64 S_pert_final = S + dS_final_requested;
            if( S_pert_final > maxS )
            {
              S_pert_final = S - dS_final_requested;                   // backward difference if at upper bound
              S_pert_final = LvArray::math::max( S_pert_final, minS );
            }
            real64 const dS_final = S_pert_final - S;
            real64 const pc_final = this->computeCapillaryPressureForSaturation(
              S, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingNonWettingCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );
            real64 const pc_pert_final = this->computeCapillaryPressureForSaturation(
              S_pert_final, mode, ipWater,
              phaseMinHistoricalVolFraction[ipWater],
              phaseMaxHistoricalVolFraction[ipWater],
              phaseMode2PeakVolFraction[ipWater],
              m_wettingNonWettingCapillaryPressureKernelWrappers,
              m_wettingCurve,
              m_nonWettingCurve,
              m_landParam[ipWater],
              m_phaseIntermediateMinVolFraction,
              m_killoughCurvatureParamCapPres,
              m_jerauldParam_a,
              m_jerauldParam_b,
              true,               // isWettingPhase
              precomputedScrt_water,
              precomputedDenomF_water,
              precomputedShy_water );
            real64 dpc_dS_final = 0.0;
            if( LvArray::math::abs( dS_final ) > 1e-14 )
            {
              dpc_dS_final = (pc_pert_final - pc_final) / dS_final;
            }
            // computeInv must return dS/dPc (not dPc/dS), matching the inverse table convention.
            // The local solver uses this as: dV/dPc = dV/dS * dS/dPc.
            real64 const dS_dPc_final = (LvArray::math::abs( dpc_dS_final ) > 1e-14) ? (1.0 / dpc_dS_final) : 0.0;
            dPhaseCapPressure_dPhaseVolFrac[ipWater][ipWater] = dS_dPc_final;
            if constexpr (ENABLE_TWOPHASE_SCANNINGINV_DEBUG) {
              std::cout << "[2P_SCAN_INV] dPc/dS=" << dpc_dS_final
                        << ", dS/dPc=" << dS_dPc_final
                        << ", Pc(S)=" << pc_final << ", Pc(S+dS)=" << pc_pert_final << std::endl;
            }
            // --- Concise one-line debug for scanning curve computeInv ---
            constexpr bool ENABLE_SCANINV_SUMMARY = false;
            if constexpr (ENABLE_SCANINV_SUMMARY) {
              std::cout << "[SCAN_INV] mode=" << static_cast< int >(mode)
                        << " Pc_tgt=" << phaseCapPressure[ipWater]
                        << " Shy=" << precomputedShy_water
                        << " Swma=" << Swma_debug
                        << " S=" << S
                        << " dS/dPc=" << dS_dPc_final
                        << (skipBisection ? " SKIP_BISECT" : " BISECT")
                        << std::endl;
            }
          }
        }                 // end if (!degenerateBracket)
      }               // End else block for scanning curves
    }
  }
}

/// kernel creation

TableCapillaryPressureHysteresis::KernelWrapper TableCapillaryPressureHysteresis::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  // Validate that the arrays are properly populated
  integer const numPhases = m_phaseNames.size();
  if( numPhases == 2 )
  {
    GEOS_THROW_IF( m_wettingNonWettingCapillaryPressureKernelWrappers.size() != 2,
                   GEOS_FMT( "{}: Expected 2 kernel wrappers for two-phase flow, but got {}. "
                             "This usually means createAllTableKernelWrappers() failed to populate the arrays. "
                             "Check that table functions '{}' and '{}' exist and are properly defined.",
                             getFullName(),
                             m_wettingNonWettingCapillaryPressureKernelWrappers.size(),
                             m_drainageWettingNonWettingCapPresTableName,
                             m_imbibitionWettingNonWettingCapPresTableName.empty() ? m_drainageWettingNonWettingCapPresTableName : m_imbibitionWettingNonWettingCapPresTableName ),
                   InputError );
  }
  else if( numPhases == 3 )
  {
    GEOS_THROW_IF( m_wettingIntermediateCapillaryPressureKernelWrappers.size() != 2,
                   GEOS_FMT( "{}: Expected 2 wetting-intermediate kernel wrappers for three-phase flow, but got {}",
                             getFullName(),
                             m_wettingIntermediateCapillaryPressureKernelWrappers.size()),
                   InputError );
    GEOS_THROW_IF( m_nonWettingIntermediateCapillaryPressureKernelWrappers.size() != 2,
                   GEOS_FMT( "{}: Expected 2 non-wetting-intermediate kernel wrappers for three-phase flow, but got {}",
                             getFullName(),
                             m_nonWettingIntermediateCapillaryPressureKernelWrappers.size()),
                   InputError );
  }

  // Validate that m_phaseHasHysteresis is properly initialized
  GEOS_THROW_IF( m_phaseHasHysteresis.size() != 2,
                 GEOS_FMT( "{}: m_phaseHasHysteresis must have size 2, but got {}",
                           getFullName(),
                           m_phaseHasHysteresis.size()),
                 InputError );

  // Validate that the historical volume fraction arrays have been resized
  // These arrays should be resized by resizeFields() before the KernelWrapper is used
  // If they're empty, it means resizeFields() hasn't been called yet
  GEOS_THROW_IF( m_phaseMaxHistoricalVolFraction.size( 0 ) == 0 || m_phaseMaxHistoricalVolFraction.size( 1 ) == 0,
                 GEOS_FMT( "{}: phaseMaxHistoricalVolFraction array has not been resized (size=[{}, {}]). "
                           "This usually means resizeFields() has not been called yet. "
                           "The arrays must be resized before the KernelWrapper can be used.",
                           getFullName(),
                           m_phaseMaxHistoricalVolFraction.size( 0 ),
                           m_phaseMaxHistoricalVolFraction.size( 1 )),
                 InputError );
  GEOS_THROW_IF( m_phaseMinHistoricalVolFraction.size( 0 ) == 0 || m_phaseMinHistoricalVolFraction.size( 1 ) == 0,
                 GEOS_FMT( "{}: phaseMinHistoricalVolFraction array has not been resized (size=[{}, {}]). "
                           "This usually means resizeFields() has not been called yet. "
                           "The arrays must be resized before the KernelWrapper can be used.",
                           getFullName(),
                           m_phaseMinHistoricalVolFraction.size( 0 ),
                           m_phaseMinHistoricalVolFraction.size( 1 )),
                 InputError );
  GEOS_THROW_IF( m_phaseMode2PeakVolFraction.size( 0 ) == 0 || m_phaseMode2PeakVolFraction.size( 1 ) == 0,
                 GEOS_FMT( "{}: phaseMode2PeakVolFraction array has not been resized (size=[{}, {}]). "
                           "This usually means resizeFields() has not been called yet. "
                           "The arrays must be resized before the KernelWrapper can be used.",
                           getFullName(),
                           m_phaseMode2PeakVolFraction.size( 0 ),
                           m_phaseMode2PeakVolFraction.size( 1 )),
                 InputError );

  // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
  return KernelWrapper( m_wettingNonWettingCapillaryPressureKernelWrappers,
                        m_inverseWettingNonWettingCapillaryPressureKernelWrappers,
                        m_wettingIntermediateCapillaryPressureKernelWrappers,
                        m_inverseWettingIntermediateCapillaryPressureKernelWrappers,
                        m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                        m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_jerauldParam_a,
                        m_jerauldParam_b,
                        m_killoughCurvatureParamCapPres,
                        m_phaseIntermediateMinVolFraction,
                        m_wettingCurve,
                        m_nonWettingCurve,
                        m_phaseMinHistoricalVolFraction,
                        m_phaseMaxHistoricalVolFraction,
                        m_phaseMode2PeakVolFraction,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_mode,
                        m_phaseTrappedVolFrac,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

void TableCapillaryPressureHysteresis::createAllTableKernelWrappers()
{
  using TPP = ThreePhasePairPhaseType;

  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();

  // Ensure m_phaseHasHysteresis is initialized before accessing it
  // This can happen if createKernelWrapper is called before postProcessInput
  if( m_phaseHasHysteresis.size() == 0 )
  {
    m_phaseHasHysteresis.resize( 2 );
    if( numPhases == 2 )
    {
      m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] =
        (m_imbibitionWettingNonWettingCapPresTableName.empty() ||
         m_imbibitionWettingNonWettingCapPresTableName == m_drainageWettingNonWettingCapPresTableName)
                        ? 0 : 1;
      m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING];
    }
    else if( numPhases == 3 )
    {
      m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] =
        (m_imbibitionWettingIntermediateCapPresTableName.empty() ||
         m_imbibitionWettingIntermediateCapPresTableName == m_drainageWettingIntermediateCapPresTableName)
                        ? 0 : 1;
      m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] =
        (m_imbibitionNonWettingIntermediateCapPresTableName.empty() ||
         m_imbibitionNonWettingIntermediateCapPresTableName == m_drainageNonWettingIntermediateCapPresTableName)
                        ? 0 : 1;
    }
  }

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_wettingNonWettingCapillaryPressureKernelWrappers.clear();
  m_inverseWettingNonWettingCapillaryPressureKernelWrappers.clear();
  m_wettingIntermediateCapillaryPressureKernelWrappers.clear();
  m_inverseWettingIntermediateCapillaryPressureKernelWrappers.clear();
  m_nonWettingIntermediateCapillaryPressureKernelWrappers.clear();
  m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers.clear();
  m_inverseTables.clear();

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( m_drainageWettingNonWettingCapPresTableName.empty(),
                   GEOS_FMT( "{}: drainageWettingNonWettingCapPressureTableName is empty for two-phase flow",
                             getFullName()),
                   InputError );

    GEOS_THROW_IF( !functionManager.hasGroup( m_drainageWettingNonWettingCapPresTableName ),
                   GEOS_FMT( "{}: the table function named '{}' could not be found",
                             getFullName(),
                             m_drainageWettingNonWettingCapPresTableName ),
                   InputError );
    TableFunction const & drainageCapPresTable = functionManager.getGroup< TableFunction >(
      m_drainageWettingNonWettingCapPresTableName );
    m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
      drainageCapPresTable.createKernelWrapper());


    string const & imbibitionTableName = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                                     ? m_imbibitionWettingNonWettingCapPresTableName
                                                                     : m_drainageWettingNonWettingCapPresTableName;
    GEOS_THROW_IF( imbibitionTableName.empty(),
                   GEOS_FMT( "{}: imbibition table name is empty for two-phase flow",
                             getFullName()),
                   InputError );
    GEOS_THROW_IF( !functionManager.hasGroup( imbibitionTableName ),
                   GEOS_FMT( "{}: the table function named '{}' could not be found",
                             getFullName(),
                             imbibitionTableName ),
                   InputError );
    TableFunction const & imbibitionWettingCapPresTable = functionManager.getGroup< TableFunction >(
      imbibitionTableName );
    m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
      imbibitionWettingCapPresTable.createKernelWrapper());

    GEOS_THROW_IF( m_wettingNonWettingCapillaryPressureKernelWrappers.size() != 2,
                   GEOS_FMT( "{}: Expected 2 kernel wrappers after creation, but got {}",
                             getFullName(),
                             m_wettingNonWettingCapillaryPressureKernelWrappers.size()),
                   InputError );

    // Create inverse tables for drainage and imbibition (for direct lookup in computeInv)
    // Inverse for drainage (index 0)
    auto const & drainageSatArrayView = drainageCapPresTable.getCoordinates()[0];
    auto const & drainagePcArrayView = drainageCapPresTable.getValues();
    std::vector< real64 > drainageSatVec( drainageSatArrayView.size());
    std::vector< real64 > drainagePcVec( drainagePcArrayView.size());
    std::copy( drainageSatArrayView.begin(), drainageSatArrayView.end(), drainageSatVec.begin());
    std::copy( drainagePcArrayView.begin(), drainagePcArrayView.end(), drainagePcVec.begin());
    // Reverse both arrays (if original Pc is decreasing in S)
    std::reverse( drainagePcVec.begin(), drainagePcVec.end());
    std::reverse( drainageSatVec.begin(), drainageSatVec.end());


    auto inverseDrainageTable = std::make_shared< TableFunction >( "inverseDrainageCapPres", this );
    real64_array invDrainagePcVec( drainagePcVec.size());
    real64_array invDrainageSatVec( drainageSatVec.size());
    std::copy( drainagePcVec.begin(), drainagePcVec.end(), invDrainagePcVec.data());
    std::copy( drainageSatVec.begin(), drainageSatVec.end(), invDrainageSatVec.data());
    array1d< real64_array > drainageCoordinates;
    drainageCoordinates.emplace_back( std::move( invDrainagePcVec ));
    std::vector< units::Unit > dimUnits = {units::Unknown};
    inverseDrainageTable->setTableCoordinates( drainageCoordinates, dimUnits );
    inverseDrainageTable->setTableValues( std::move( invDrainageSatVec ), units::Unknown );
    inverseDrainageTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseWettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
      inverseDrainageTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseDrainageTable ));

    // Inverse for imbibition (index 1)
    auto const & imbibitionSatArrayView = imbibitionWettingCapPresTable.getCoordinates()[0];
    auto const & imbibitionPcArrayView = imbibitionWettingCapPresTable.getValues();
    std::vector< real64 > imbibitionSatVec( imbibitionSatArrayView.size());
    std::vector< real64 > imbibitionPcVec( imbibitionPcArrayView.size());
    std::copy( imbibitionSatArrayView.begin(), imbibitionSatArrayView.end(), imbibitionSatVec.begin());
    std::copy( imbibitionPcArrayView.begin(), imbibitionPcArrayView.end(), imbibitionPcVec.begin());
    // Reverse both arrays (if original Pc is decreasing in S)
    std::reverse( imbibitionPcVec.begin(), imbibitionPcVec.end());
    std::reverse( imbibitionSatVec.begin(), imbibitionSatVec.end());

    auto inverseImbibitionTable = std::make_shared< TableFunction >( "inverseImbibitionCapPres", this );
    real64_array invImbibitionPcVec( imbibitionPcVec.size());
    real64_array invImbibitionSatVec( imbibitionSatVec.size());
    std::copy( imbibitionPcVec.begin(), imbibitionPcVec.end(), invImbibitionPcVec.data());
    std::copy( imbibitionSatVec.begin(), imbibitionSatVec.end(), invImbibitionSatVec.data());
    array1d< real64_array > imbibitionCoordinates;
    imbibitionCoordinates.emplace_back( std::move( invImbibitionPcVec ));
    inverseImbibitionTable->setTableCoordinates( imbibitionCoordinates, dimUnits );
    inverseImbibitionTable->setTableValues( std::move( invImbibitionSatVec ), units::Unknown );
    inverseImbibitionTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseWettingNonWettingCapillaryPressureKernelWrappers.emplace_back(
      inverseImbibitionTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseImbibitionTable ));

  }
  else if( numPhases == 3 )
  {
    TableFunction const & drainageWICapPres = functionManager.getGroup< TableFunction >(
      m_drainageWettingIntermediateCapPresTableName );
    m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      drainageWICapPres.createKernelWrapper());

    TableFunction const & drainageNWICapPres = functionManager.getGroup< TableFunction >(
      m_drainageNonWettingIntermediateCapPresTableName );
    m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      drainageNWICapPres.createKernelWrapper());

    TableFunction const & imbibitionWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                           ? functionManager.getGroup< TableFunction >(
      m_imbibitionWettingIntermediateCapPresTableName )
                                                           : functionManager.getGroup< TableFunction >(
      m_drainageWettingIntermediateCapPresTableName );
    m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      imbibitionWICapPres.createKernelWrapper());

    TableFunction const & imbibitionNWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING]
                                                            ? functionManager.getGroup< TableFunction >(
      m_imbibitionNonWettingIntermediateCapPresTableName )
                                                            : functionManager.getGroup< TableFunction >(
      m_drainageNonWettingIntermediateCapPresTableName );
    m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      imbibitionNWICapPres.createKernelWrapper());

    // Create inverse tables for 3-phase (wetting-intermediate and non-wetting-intermediate)
    // Inverse for wetting-intermediate drainage (index 0)
    auto const & wiDrainageSatArrayView = drainageWICapPres.getCoordinates()[0];
    auto const & wiDrainagePcArrayView = drainageWICapPres.getValues();
    std::vector< real64 > wiDrainageSatVec( wiDrainageSatArrayView.size());
    std::vector< real64 > wiDrainagePcVec( wiDrainagePcArrayView.size());
    std::copy( wiDrainageSatArrayView.begin(), wiDrainageSatArrayView.end(), wiDrainageSatVec.begin());
    std::copy( wiDrainagePcArrayView.begin(), wiDrainagePcArrayView.end(), wiDrainagePcVec.begin());
    std::reverse( wiDrainagePcVec.begin(), wiDrainagePcVec.end());
    std::reverse( wiDrainageSatVec.begin(), wiDrainageSatVec.end());

    auto inverseWIDrainageTable = std::make_shared< TableFunction >( "inverseWIDrainageCapPres", this );
    real64_array invWIDrainagePcVec( wiDrainagePcVec.size());
    real64_array invWIDrainageSatVec( wiDrainageSatVec.size());
    std::copy( wiDrainagePcVec.begin(), wiDrainagePcVec.end(), invWIDrainagePcVec.data());
    std::copy( wiDrainageSatVec.begin(), wiDrainageSatVec.end(), invWIDrainageSatVec.data());
    array1d< real64_array > wiDrainageCoordinates;
    wiDrainageCoordinates.emplace_back( std::move( invWIDrainagePcVec ));
    std::vector< units::Unit > dimUnits = {units::Unknown};
    inverseWIDrainageTable->setTableCoordinates( wiDrainageCoordinates, dimUnits );
    inverseWIDrainageTable->setTableValues( std::move( invWIDrainageSatVec ), units::Unknown );
    inverseWIDrainageTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      inverseWIDrainageTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseWIDrainageTable ));

    // Inverse for wetting-intermediate imbibition (index 1)
    auto const & wiImbibitionSatArrayView = imbibitionWICapPres.getCoordinates()[0];
    auto const & wiImbibitionPcArrayView = imbibitionWICapPres.getValues();
    std::vector< real64 > wiImbibitionSatVec( wiImbibitionSatArrayView.size());
    std::vector< real64 > wiImbibitionPcVec( wiImbibitionPcArrayView.size());
    std::copy( wiImbibitionSatArrayView.begin(), wiImbibitionSatArrayView.end(), wiImbibitionSatVec.begin());
    std::copy( wiImbibitionPcArrayView.begin(), wiImbibitionPcArrayView.end(), wiImbibitionPcVec.begin());
    std::reverse( wiImbibitionPcVec.begin(), wiImbibitionPcVec.end());
    std::reverse( wiImbibitionSatVec.begin(), wiImbibitionSatVec.end());

    auto inverseWIImbibitionTable = std::make_shared< TableFunction >( "inverseWIImbibitionCapPres", this );
    real64_array invWIImbibitionPcVec( wiImbibitionPcVec.size());
    real64_array invWIImbibitionSatVec( wiImbibitionSatVec.size());
    std::copy( wiImbibitionPcVec.begin(), wiImbibitionPcVec.end(), invWIImbibitionPcVec.data());
    std::copy( wiImbibitionSatVec.begin(), wiImbibitionSatVec.end(), invWIImbibitionSatVec.data());
    array1d< real64_array > wiImbibitionCoordinates;
    wiImbibitionCoordinates.emplace_back( std::move( invWIImbibitionPcVec ));
    inverseWIImbibitionTable->setTableCoordinates( wiImbibitionCoordinates, dimUnits );
    inverseWIImbibitionTable->setTableValues( std::move( invWIImbibitionSatVec ), units::Unknown );
    inverseWIImbibitionTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      inverseWIImbibitionTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseWIImbibitionTable ));

    // Inverse for non-wetting-intermediate drainage (index 0)
    auto const & nwiDrainageSatArrayView = drainageNWICapPres.getCoordinates()[0];
    auto const & nwiDrainagePcArrayView = drainageNWICapPres.getValues();
    std::vector< real64 > nwiDrainageSatVec( nwiDrainageSatArrayView.size());
    std::vector< real64 > nwiDrainagePcVec( nwiDrainagePcArrayView.size());
    std::copy( nwiDrainageSatArrayView.begin(), nwiDrainageSatArrayView.end(), nwiDrainageSatVec.begin());
    std::copy( nwiDrainagePcArrayView.begin(), nwiDrainagePcArrayView.end(), nwiDrainagePcVec.begin());
    std::reverse( nwiDrainagePcVec.begin(), nwiDrainagePcVec.end());
    std::reverse( nwiDrainageSatVec.begin(), nwiDrainageSatVec.end());

    auto inverseNWIDrainageTable = std::make_shared< TableFunction >( "inverseNWIDrainageCapPres", this );
    real64_array invNWIDrainagePcVec( nwiDrainagePcVec.size());
    real64_array invNWIDrainageSatVec( nwiDrainageSatVec.size());
    std::copy( nwiDrainagePcVec.begin(), nwiDrainagePcVec.end(), invNWIDrainagePcVec.data());
    std::copy( nwiDrainageSatVec.begin(), nwiDrainageSatVec.end(), invNWIDrainageSatVec.data());
    array1d< real64_array > nwiDrainageCoordinates;
    nwiDrainageCoordinates.emplace_back( std::move( invNWIDrainagePcVec ));
    inverseNWIDrainageTable->setTableCoordinates( nwiDrainageCoordinates, dimUnits );
    inverseNWIDrainageTable->setTableValues( std::move( invNWIDrainageSatVec ), units::Unknown );
    inverseNWIDrainageTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      inverseNWIDrainageTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseNWIDrainageTable ));

    // Inverse for non-wetting-intermediate imbibition (index 1)
    auto const & nwiImbibitionSatArrayView = imbibitionNWICapPres.getCoordinates()[0];
    auto const & nwiImbibitionPcArrayView = imbibitionNWICapPres.getValues();
    std::vector< real64 > nwiImbibitionSatVec( nwiImbibitionSatArrayView.size());
    std::vector< real64 > nwiImbibitionPcVec( nwiImbibitionPcArrayView.size());
    std::copy( nwiImbibitionSatArrayView.begin(), nwiImbibitionSatArrayView.end(), nwiImbibitionSatVec.begin());
    std::copy( nwiImbibitionPcArrayView.begin(), nwiImbibitionPcArrayView.end(), nwiImbibitionPcVec.begin());
    std::reverse( nwiImbibitionPcVec.begin(), nwiImbibitionPcVec.end());
    std::reverse( nwiImbibitionSatVec.begin(), nwiImbibitionSatVec.end());

    auto inverseNWIImbibitionTable = std::make_shared< TableFunction >( "inverseNWIImbibitionCapPres", this );
    real64_array invNWIImbibitionPcVec( nwiImbibitionPcVec.size());
    real64_array invNWIImbibitionSatVec( nwiImbibitionSatVec.size());
    std::copy( nwiImbibitionPcVec.begin(), nwiImbibitionPcVec.end(), invNWIImbibitionPcVec.data());
    std::copy( nwiImbibitionSatVec.begin(), nwiImbibitionSatVec.end(), invNWIImbibitionSatVec.data());
    array1d< real64_array > nwiImbibitionCoordinates;
    nwiImbibitionCoordinates.emplace_back( std::move( invNWIImbibitionPcVec ));
    inverseNWIImbibitionTable->setTableCoordinates( nwiImbibitionCoordinates, dimUnits );
    inverseNWIImbibitionTable->setTableValues( std::move( invNWIImbibitionSatVec ), units::Unknown );
    inverseNWIImbibitionTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back(
      inverseNWIImbibitionTable->createKernelWrapper());
    m_inverseTables.emplace_back( std::move( inverseNWIImbibitionTable ));
  }

}

///kernel ctor
TableCapillaryPressureHysteresis::KernelWrapper::KernelWrapper(
  arrayView1d< const TableFunction::KernelWrapper > const & wettingNonWettingCapillaryPressureKernelWrappers,
  arrayView1d< const TableFunction::KernelWrapper > const & inverseWettingNonWettingCapillaryPressureKernelWrappers,
  arrayView1d< const TableFunction::KernelWrapper > const & wettingIntermediateCapillaryPressureKernelWrappers,
  arrayView1d< const TableFunction::KernelWrapper > const & inverseWettingIntermediateCapillaryPressureKernelWrappers,
  arrayView1d< const TableFunction::KernelWrapper > const & nonWettingIntermediateCapillaryPressureKernelWrappers,
  arrayView1d< const TableFunction::KernelWrapper > const & inverseNonWettingIntermediateCapillaryPressureKernelWrappers,
  const arrayView1d< const geos::integer > & phaseHasHysteresis,
  const arrayView1d< const geos::real64 > & landParam,
  const real64 & jerauldParam_a,
  const real64 & jerauldParam_b,
  const real64 & killoughCurvaturePcParam,
  const geos::real64 & phaseIntermediateMinVolFraction,
  const KilloughHysteresis::HysteresisCurve & wettingCurve,
  const KilloughHysteresis::HysteresisCurve & nonWettingCurve,
  const arrayView2d< const geos::real64, compflow::USD_PHASE > & phaseMinHistoricalVolFraction,
  const arrayView2d< const geos::real64, compflow::USD_PHASE > & phaseMaxHistoricalVolFraction,
  arrayView2d< geos::real64, compflow::USD_PHASE > & phaseMode2PeakVolFraction,
  arrayView1d< integer const > const & phaseTypes,
  arrayView1d< integer const > const & phaseOrder,
  arrayView1d< integer > const & mode,
  arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrapped,
  const arrayView3d< geos::real64, relperm::USD_RELPERM > & phaseCapPressure,
  const arrayView4d< geos::real64, relperm::USD_RELPERM_DS > & dPhaseCapPressure_dPhaseVolFrac )
  :
  CapillaryPressureBaseUpdate( phaseTypes,
                               phaseOrder,
                               phaseTrapped,
                               phaseCapPressure,
                               dPhaseCapPressure_dPhaseVolFrac ),
  m_wettingNonWettingCapillaryPressureKernelWrappers( wettingNonWettingCapillaryPressureKernelWrappers ),
  m_inverseWettingNonWettingCapillaryPressureKernelWrappers( inverseWettingNonWettingCapillaryPressureKernelWrappers ),
  m_wettingIntermediateCapillaryPressureKernelWrappers(
    wettingIntermediateCapillaryPressureKernelWrappers ),
  m_inverseWettingIntermediateCapillaryPressureKernelWrappers( inverseWettingIntermediateCapillaryPressureKernelWrappers ),
  m_nonWettingIntermediateCapillaryPressureKernelWrappers(
    nonWettingIntermediateCapillaryPressureKernelWrappers ),
  m_inverseNonWettingIntermediateCapillaryPressureKernelWrappers( inverseNonWettingIntermediateCapillaryPressureKernelWrappers ),
  m_phaseHasHysteresis( phaseHasHysteresis ),
  m_landParam( landParam ),
  m_jerauldParam_a( jerauldParam_a ),
  m_jerauldParam_b( jerauldParam_b ),
  m_killoughCurvatureParamCapPres( killoughCurvaturePcParam ),
  m_phaseIntermediateMinVolFraction( phaseIntermediateMinVolFraction ),
  m_wettingCurve( wettingCurve ),
  m_nonWettingCurve( nonWettingCurve ),
  m_phaseMinHistoricalVolFraction( phaseMinHistoricalVolFraction ),
  m_phaseMaxHistoricalVolFraction( phaseMaxHistoricalVolFraction ),
  m_phaseMode2PeakVolFraction( phaseMode2PeakVolFraction ),
  m_mode( mode ) {}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableCapillaryPressureHysteresis, std::string const &, Group * const )

}     // namespace constitutive
} // namespace geos
