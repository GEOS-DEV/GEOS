/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SetInitialTemperatureAndPressureMPMEvent.hpp
 */

#ifndef GEOSX_SETINITIALTEMPERATUREANDPRESSURE_MPMEVENT_HPP_
#define GEOSX_SETINITIALTEMPERATUREANDPRESSURE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class SetInitialTemperatureAndPressureMPMEvent
 *
 * Initialization event that asks each selected material model to compute an EOS-consistent
 * hydrostatic state at a user-specified pressure and temperature, then lets the MPM solver
 * copy the returned density, energy, stress, and kinematic reference state to particles.
 */
class SetInitialTemperatureAndPressureMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  SetInitialTemperatureAndPressureMPMEvent( const string & name,
                                            Group * const parent );

  /// Destructor
  virtual ~SetInitialTemperatureAndPressureMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "SetInitialTemperatureAndPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * pressureString() { return "pressure"; }
    static constexpr char const * temperatureString() { return "temperature"; }
    static constexpr char const * targetRegionsString() { return "targetRegions"; }
    static constexpr char const * initialPhaseString() { return "initialPhase"; }
    static constexpr char const * initializeStressString() { return "initializeStress"; }
  } SetInitialTemperatureAndPressureMPMEventViewKeys;
  /// @endcond

  real64 getPressure() const { return m_pressure; }
  real64 getTemperature() const { return m_temperature; }
  string_array const & getTargetRegions() const { return m_targetRegions; }
  string const & getInitialPhase() const { return m_initialPhase; }
  integer getInitializeStress() const { return m_initializeStress; }

protected:
  virtual void postInputInitialization() override final;

  real64 m_pressure;
  real64 m_temperature;
  string_array m_targetRegions;
  string m_initialPhase;
  integer m_initializeStress;
};

} /* namespace geos */

#endif /* GEOSX_SETINITIALTEMPERATUREANDPRESSURE_MPMEVENT_HPP_ */
