/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HaltEvent.hpp
 */


#ifndef GEOSX_MANAGERS_EVENTS_HALTEVENT_HPP_
#define GEOSX_MANAGERS_EVENTS_HALTEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

/**
 * @class HaltEvent
 * An event type that is designed to look at the external clock.
 * This is useful for managing wall time limitations.
 */
class HaltEvent : public EventBase
{
public:
  /// Constructor
  HaltEvent( const std::string & name,
             Group * const parent );

  /// Destructor
  virtual ~HaltEvent() override;

  // Catalog name interface
  static string CatalogName() { return "HaltEvent"; }

  /**
   * This event is designed to look at the external clock. Currently,
   * if the event is triggered it will set a flag, which will
   * instruct the code to exit.  This is useful for managing walltime
   */
  virtual void EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    dataRepository::Group * domain ) override;

  /// Timing values
  real64 m_externalStartTime;
  real64 m_externalLastTime;
  real64 m_externalDt;
  real64 m_maxRuntime;

  struct viewKeyStruct
  {
    static constexpr auto maxRuntimeString = "maxRuntime";

    dataRepository::ViewKey maxRuntime = { "maxRuntime" };
  } haltEventViewKeys;

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_HALTEVENT_HPP_ */
