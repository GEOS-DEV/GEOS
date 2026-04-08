/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CO2BrineFluid.cpp
 */
#include "CO2BrineFluid.hpp"

#include "constitutive/fluid/multifluid/LogLevelsInfo.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "common/Units.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename PHASE1, typename PHASE2, typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PhillipsBrineDensity, PhillipsBrineViscosity, NoOpPVTFunction >,
                        PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, NoOpPVTFunction >,
                        CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< PhillipsBrineDensity, PhillipsBrineViscosity, BrineEnthalpy >,
                        PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, CO2Enthalpy >,
                        CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsThermalFluid"; }
};

template<> class
  TwoPhaseCatalogNames< PhaseModel< EzrokhiBrineDensity, EzrokhiBrineViscosity, NoOpPVTFunction >,
                        PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, NoOpPVTFunction >,
                        CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< EzrokhiBrineDensity, EzrokhiBrineViscosity, BrineEnthalpy >,
                        PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, CO2Enthalpy >,
                        CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiThermalFluid"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename PHASE1, typename PHASE2, typename FLASH >
string CO2BrineFluid< PHASE1, PHASE2, FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< PHASE1, PHASE2, FLASH >::name();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::
CO2BrineFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "When set to 1, write PVT tables into a CSV file.\n"
                    "if the table is requested to be output in the log, and it is too large, a CSV file will be generated even if `writeCSV` is set to 0." ).
    setDefaultValue( 0 );

  registerWrapper( viewKeyStruct::checkPhasePresenceString(), &m_checkPhasePresence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Check phase presence when computing density and viscosity" ).
    setApplyDefaultValue( 0 );

  // Attach the fluid properties
  bool constexpr EZROKHI_DENSITY = std::is_same_v< typename PHASE1::Density, EzrokhiBrineDensity >;
  bool constexpr EZROKHI_VISCOSITY = std::is_same_v< typename PHASE1::Viscosity, EzrokhiBrineViscosity >;
  m_brineFluidParameters.registerOnFluid< true, EZROKHI_DENSITY, EZROKHI_VISCOSITY >( this );

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if constexpr ( isThermalType() )
  {
    getField< fields::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getField< fields::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }

  addLogLevel< logInfo::TableLogOutput >();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
std::unique_ptr< ConstitutiveBase >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  CO2BrineFluid & newConstitutiveRelation = dynamicCast< CO2BrineFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename PHASE1, typename PHASE2, typename FLASH >
integer CO2BrineFluid< PHASE1, PHASE2, FLASH >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::checkTablesParameters( real64 const pressure,
                                                                    real64 const temperature ) const
{
  if( !m_checkPVTTablesRanges )
  {
    return;
  }

  real64 const temperatureInCelsius = units::convertKToC( temperature );
  try
  {
    m_phase1->density.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase1->viscosity.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase1->enthalpy.checkTablesParameters( pressure, temperatureInCelsius );
  }
  catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error for {} phase.", m_phaseNames[m_p1Index] );
    ErrorLogger::global().modifyCurrentExceptionMessage()
      .addToMsg( errorMsg )
      .addContextInfo( getDataContext().getContextInfo().setPriority( 2 ) );
    throw SimulationError( ex, errorMsg );
  }

  try
  {
    m_phase2->density.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase2->viscosity.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase2->enthalpy.checkTablesParameters( pressure, temperatureInCelsius );
  }
  catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error for {} phase.", m_phaseNames[m_p2Index] );
    ErrorLogger::global().modifyCurrentExceptionMessage()
      .addToMsg( errorMsg )
      .addContextInfo( getDataContext().getContextInfo().setPriority( 2 ) );
    throw SimulationError( ex, errorMsg );
  }

  try
  {
    m_flash->checkTablesParameters( pressure, temperatureInCelsius );
  }
  catch( SimulationError const & ex )
  {
    string const errorMsg = "Table input error for flash parameters";
    ErrorLogger::global().modifyCurrentExceptionMessage()
      .addToMsg( errorMsg )
      .addContextInfo( getDataContext().getContextInfo().setPriority( 2 ) );
    throw SimulationError( ex, errorMsg );
  }
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::initializePreSubGroups()
{
#if defined(GEOS_DEVICE_COMPILE)
  if constexpr (std::is_same_v< CO2BrineFluid< PHASE1, PHASE2, FLASH >, CO2BrineEzrokhiThermalFluid >)
  {
    GEOS_THROW( GEOS_FMT( "The `{}` model is disabled for now. Please use the other thermal CO2-brine model instead: `{}`",
                          CO2BrineEzrokhiThermalFluid::catalogName(),
                          CO2BrinePhillipsThermalFluid::catalogName() ),
                InputError, getDataContext() );
  }
#endif
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( numFluidPhases(), 2,
                        "invalid number of phases",
                        InputError, getDataContext() );
  GEOS_THROW_IF_NE_MSG( numFluidComponents(), 2,
                        "invalid number of components",
                        InputError, getDataContext() );

  // NOTE: for now, the names of the phases are still hardcoded here
  // Later, we could read them from the XML file and we would then have a general class here

  string const expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_p1Index = PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );

  string const expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_p2Index = PVTFunctionHelpers::findName( m_phaseNames, expectedGasPhaseNames, viewKeyStruct::phaseNamesString() );

  // Validate the brine fluid properties
  bool constexpr EZROKHI_DENSITY = std::is_same_v< typename PHASE1::Density, EzrokhiBrineDensity >;
  bool constexpr EZROKHI_VISCOSITY = std::is_same_v< typename PHASE1::Viscosity, EzrokhiBrineViscosity >;
  m_brineFluidParameters.postInputInitialization< true, EZROKHI_DENSITY, EZROKHI_VISCOSITY >( this );

  createPVTModels();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::createPVTModels()
{
  // then, we are ready to instantiate the phase models
  bool const isClone = this->isClone();
  TableFunction::OutputOptions const outputOpts = {
    !isClone && m_writeCSV, // writeCSV
    !isClone && isLogLevelActive< logInfo::TableLogOutput >( this->getLogLevel()) // writeInLog
  };

  // 1) Create phase PVT models
  m_phase1 = std::make_unique< PHASE1 >( getName() + "_phaseModel1",
                                         m_brineFluidParameters,
                                         m_componentNames,
                                         m_componentMolarWeight,
                                         outputOpts );
  m_phase2 = std::make_unique< PHASE2 >( getName() + "_phaseModel2",
                                         m_brineFluidParameters,
                                         m_componentNames,
                                         m_componentMolarWeight,
                                         outputOpts );

  // 2) Create the flash model
  m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                       m_brineFluidParameters,
                                       m_phaseNames,
                                       m_componentNames,
                                       m_componentMolarWeight,
                                       outputOpts );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
typename CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper
CO2BrineFluid< PHASE1, PHASE2, FLASH >::createKernelWrapper()
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_phase1,
                        *m_phase2,
                        *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_checkPhasePresence,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper::
  KernelWrapper( integer const p1Index,
                 integer const p2Index,
                 PHASE1 const & phase1,
                 PHASE2 const & phase2,
                 FLASH const & flash,
                 arrayView1d< geos::real64 const > componentMolarWeight,
                 bool const useMass,
                 bool const isThermal,
                 bool const checkPhasePresence,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_p1Index( p1Index ),
  m_p2Index( p2Index ),
  m_isThermal( isThermal ),
  m_checkPhasePresence( checkPhasePresence ),
  m_phase1( phase1.createKernelWrapper() ),
  m_phase2( phase2.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class CO2BrineFluid< PhaseModel< PhillipsBrineDensity, PhillipsBrineViscosity, NoOpPVTFunction >,
                              PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, NoOpPVTFunction >,
                              CO2Solubility >;
template class CO2BrineFluid< PhaseModel< PhillipsBrineDensity, PhillipsBrineViscosity, BrineEnthalpy >,
                              PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, CO2Enthalpy >,
                              CO2Solubility >;

template class CO2BrineFluid< PhaseModel< EzrokhiBrineDensity, EzrokhiBrineViscosity, NoOpPVTFunction >,
                              PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, NoOpPVTFunction >,
                              CO2Solubility >;
template class CO2BrineFluid< PhaseModel< EzrokhiBrineDensity, EzrokhiBrineViscosity, BrineEnthalpy >,
                              PhaseModel< SpanWagnerCO2Density, FenghourCO2Viscosity, CO2Enthalpy >,
                              CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsThermalFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiThermalFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geos
