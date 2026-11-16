/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_THERMAL_SETTINGS_HPP_
#define GEOS_MPM_THERMAL_SETTINGS_HPP_

#include "dataRepository/Group.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace geos
{

/** Selectors in different categories are intersected; values in one list are ORed.
 * Membership and initial mass are cached on particles and inherited by split children.
 */
class MPMThermalParticleCondition : public dataRepository::Group
{
public:
  MPMThermalParticleCondition( string const & name, Group * parent, bool const isothermal );
  bool matches( string const & region, globalIndex id, integer group, integer color ) const;
  bool active( real64 time ) const { return time >= startTime && time < endTime; }
  real64 targetTemperature( real64 time ) const
  { return temperature + temperatureRate * ( time - startTime ); }
  void postInputInitialization() override;

  bool const isothermal;
  array1d< globalIndex > particleIDs;
  array1d< integer > contactGroups;
  array1d< integer > particleColors;
  string_array targetRegions;
  real64 startTime = 0.0;
  real64 endTime = std::numeric_limits< real64 >::max();
  real64 power = 0.0;
  real64 specificPower = 0.0;
  real64 volumetricPower = 0.0;
  real64 temperature = 293.15;
  real64 temperatureRate = 0.0;
};

class MPMThermalContactRule : public dataRepository::Group
{
public:
  MPMThermalContactRule( string const & name, Group * parent );
  void postInputInitialization() override;
  bool matches( string const & a, string const & b, integer ga, integer gb ) const;
  integer priority() const { return ( contactGroups.size() > 0 ? 2 : 0 ) + ( materials.size() > 0 ? 1 : 0 ); }

  string_array materials;
  array1d< integer > contactGroups;
  string mode = "Resistance";
  real64 resistance = 0.0;
};

/** Optional <Thermal> child of SolidMechanics_MPM. All state is restartable. */
class MPMThermalSettings : public dataRepository::Group
{
public:
  MPMThermalSettings( string const & name, Group * parent );
  Group * createChild( string const & key, string const & name ) override;
  void expandObjectCatalogs() override;
  void postInputInitialization() override;
  real64 nextConditionTime( real64 time ) const;
  MPMThermalContactRule const * contactRule( string const & a, string const & b,
                                            integer ga, integer gb ) const;

  integer enabled = 0;
  integer energyFormulationVersion = 0;
  string defaultContactMode = "Resistance";
  real64 defaultContactResistance = 0.0;
  real64 gapTolerance = 1.0e-6;
  real64 safetyFactor = 0.45;
  integer maxSubsteps = 10000;
  integer areaQuadrature = 12;
  real64 minimumCapacityFraction = 1.0e-14;
  integer diagnostics = 0;
  string mechanicalCoupling = "Coupled";
  integer maxCouplingIterations = 100;
  real64 couplingTolerance = 1.0e-10;
  real64 cumulativeMechanicalWork = 0.0;
  real64 cumulativePlasticHeat = 0.0;
  real64 suggestedTimeStep = std::numeric_limits< real64 >::max();
  real64 cumulativeSourceEnergy = 0.0;
  real64 cumulativeReservoirEnergy = 0.0;
  string_array conditionNames;
  string_array materialNames;
};

} // namespace geos
#endif
