/* SPDX-License-Identifier: LGPL-2.1-only */
#include "MPMThermalSettings.hpp"

namespace geos
{
using namespace dataRepository;
namespace
{
template< typename LIST, typename T >
bool contains( LIST const & list, T const & value )
{ return std::find( list.begin(), list.end(), value ) != list.end(); }

bool validMode( string const & mode )
{ return mode == "Insulated" || mode == "Perfect" || mode == "Resistance"; }
}

MPMThermalParticleCondition::MPMThermalParticleCondition( string const & name, Group * parent, bool const fixed ):
  Group( name, parent ), isothermal( fixed )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  registerWrapper( "particleIDs", &particleIDs ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "contactGroups", &contactGroups ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "particleColors", &particleColors ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "targetRegions", &targetRegions ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "startTime", &startTime ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( startTime );
  registerWrapper( "endTime", &endTime ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( endTime );
  if( fixed )
  {
    registerWrapper( "temperature", &temperature ).setInputFlag( InputFlags::REQUIRED );
    registerWrapper( "temperatureRate", &temperatureRate ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
  }
  else
  {
    registerWrapper( "power", &power ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 ).
      setDescription( "Power per original selected particle; split children receive their mass fraction." );
    registerWrapper( "specificPower", &specificPower ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
    registerWrapper( "volumetricPower", &volumetricPower ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
  }
}

bool MPMThermalParticleCondition::matches( string const & region, globalIndex const id,
                                          integer const group, integer const color ) const
{
  return ( targetRegions.empty() || contains( targetRegions, region ) ) &&
         ( particleIDs.empty() || contains( particleIDs, id ) ) &&
         ( contactGroups.empty() || contains( contactGroups, group ) ) &&
         ( particleColors.empty() || contains( particleColors, color ) );
}

void MPMThermalParticleCondition::postInputInitialization()
{
  Group::postInputInitialization();
  GEOS_ERROR_IF( !std::isfinite( startTime ) || startTime < 0.0 ||
                 !std::isfinite( endTime ) || !( endTime > startTime ) ||
                 !std::isfinite( power ) || !std::isfinite( specificPower ) || !std::isfinite( volumetricPower ) ||
                 !std::isfinite( temperature ) || temperature < 0.0 || !std::isfinite( temperatureRate ),
                 getDataContext() << ": invalid particle thermal condition" );
}

MPMThermalContactRule::MPMThermalContactRule( string const & name, Group * parent ): Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  registerWrapper( "materials", &materials ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "contactGroups", &contactGroups ).setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "mode", &mode ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( mode );
  registerWrapper( "resistance", &resistance ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 ).
    setDescription( "Area-specific thermal contact resistance [area K/power]" );
}

void MPMThermalContactRule::postInputInitialization()
{
  Group::postInputInitialization();
  GEOS_ERROR_IF( ( materials.size() != 0 && materials.size() != 2 ) ||
                 ( contactGroups.size() != 0 && contactGroups.size() != 2 ) || priority() == 0 ||
                 !validMode( mode ) || !std::isfinite( resistance ) || resistance < 0.0,
                 getDataContext() << ": contact rules require a material pair and/or group pair and a valid law" );
  for( integer const g : contactGroups ) GEOS_ERROR_IF( g < 0, "Thermal contact groups must be nonnegative." );
}

bool MPMThermalContactRule::matches( string const & a, string const & b, integer const ga, integer const gb ) const
{
  // Keep the material/group association on each side when both selectors are supplied.
  bool const forward = ( materials.empty() || ( materials[0] == a && materials[1] == b ) ) &&
                       ( contactGroups.empty() || ( contactGroups[0] == ga && contactGroups[1] == gb ) );
  bool const reverse = ( materials.empty() || ( materials[0] == b && materials[1] == a ) ) &&
                       ( contactGroups.empty() || ( contactGroups[0] == gb && contactGroups[1] == ga ) );
  return forward || reverse;
}

MPMThermalSettings::MPMThermalSettings( string const & name, Group * parent ): Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  registerWrapper( "energyFormulationVersion", &energyFormulationVersion ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "enabled", &enabled ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0 );
  registerWrapper( "defaultContactMode", &defaultContactMode ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( defaultContactMode );
  registerWrapper( "defaultContactResistance", &defaultContactResistance ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
  registerWrapper( "gapTolerance", &gapTolerance ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( gapTolerance ).
    setDescription( "Maximum reconstructed gap for solid-contact heat transfer, as a fraction of minimum active cell size." );
  registerWrapper( "safetyFactor", &safetyFactor ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( safetyFactor );
  registerWrapper( "maxSubsteps", &maxSubsteps ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( maxSubsteps );
  registerWrapper( "areaQuadrature", &areaQuadrature ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( areaQuadrature );
  registerWrapper( "minimumCapacityFraction", &minimumCapacityFraction ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( minimumCapacityFraction );
  registerWrapper( "mechanicalCoupling", &mechanicalCoupling ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( mechanicalCoupling ).
    setDescription( "Coupled solves particle temperature and plastic return together; Lagged freezes the yield temperature for one mechanical increment." );
  registerWrapper( "maxCouplingIterations", &maxCouplingIterations ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( maxCouplingIterations );
  registerWrapper( "couplingTolerance", &couplingTolerance ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( couplingTolerance );
  registerWrapper( "cumulativeMechanicalWork", &cumulativeMechanicalWork ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "cumulativePlasticHeat", &cumulativePlasticHeat ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "diagnostics", &diagnostics ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0 );
  registerWrapper( "suggestedTimeStep", &suggestedTimeStep ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "cumulativeSourceEnergy", &cumulativeSourceEnergy ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "cumulativeReservoirEnergy", &cumulativeReservoirEnergy ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "conditionNames", &conditionNames ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
  registerWrapper( "materialNames", &materialNames ).setInputFlag( InputFlags::FALSE ).setRestartFlags( RestartFlags::WRITE_AND_READ );
}

Group * MPMThermalSettings::createChild( string const & key, string const & name )
{
  if( key == "HeatSource" || key == "Isothermal" )
    return &registerGroup( name, std::make_unique< MPMThermalParticleCondition >( name, this, key == "Isothermal" ) );
  if( key == "ContactResistance" )
    return &registerGroup( name, std::make_unique< MPMThermalContactRule >( name, this ) );
  GEOS_ERROR( "Unknown Thermal child '" << key << "'. Expected HeatSource, Isothermal, or ContactResistance." );
  return nullptr;
}

void MPMThermalSettings::expandObjectCatalogs()
{
  createChild( "HeatSource", "HeatSource" );
  createChild( "Isothermal", "Isothermal" );
  createChild( "ContactResistance", "ContactResistance" );
}

void MPMThermalSettings::postInputInitialization()
{
  Group::postInputInitialization();
  GEOS_ERROR_IF( ( enabled != 0 && enabled != 1 ) || !validMode( defaultContactMode ) ||
                 !std::isfinite( defaultContactResistance ) || defaultContactResistance < 0.0 ||
                 !std::isfinite( gapTolerance ) || gapTolerance < 0.0 ||
                 !std::isfinite( safetyFactor ) || !( safetyFactor > 0.0 && safetyFactor < 1.0 ) ||
                 !std::isfinite( minimumCapacityFraction ) || !( minimumCapacityFraction > 0.0 && minimumCapacityFraction < 1.0 ) ||
                 maxSubsteps < 1 || areaQuadrature < 2 ||
                 ( mechanicalCoupling != "Coupled" && mechanicalCoupling != "Lagged" ) ||
                 maxCouplingIterations < 2 || !std::isfinite( couplingTolerance ) || !( couplingTolerance > 0.0 && couplingTolerance < 1.0 ),
                 getDataContext() << ": invalid thermal solver controls" );
}

real64 MPMThermalSettings::nextConditionTime( real64 const time ) const
{
  real64 next = std::numeric_limits< real64 >::max();
  forSubGroups< MPMThermalParticleCondition >( [&]( MPMThermalParticleCondition const & c )
  {
    if( c.startTime > time ) next = std::min( next, c.startTime );
    if( c.endTime > time ) next = std::min( next, c.endTime );
  } );
  return next;
}

MPMThermalContactRule const * MPMThermalSettings::contactRule( string const & a, string const & b,
                                                            integer const ga, integer const gb ) const
{
  MPMThermalContactRule const * selected = nullptr;
  forSubGroups< MPMThermalContactRule >( [&]( MPMThermalContactRule const & rule )
  {
    if( !rule.matches( a, b, ga, gb ) ) return;
    GEOS_ERROR_IF( selected != nullptr && selected->priority() == rule.priority(),
                   "Ambiguous thermal contact rules '" << selected->getName() << "' and '" << rule.getName() << "'." );
    if( selected == nullptr || rule.priority() > selected->priority() ) selected = &rule;
  } );
  return selected;
}

} // namespace geos
