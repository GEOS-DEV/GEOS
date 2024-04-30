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

#ifndef GEOS_RESIDUALDUMPEVENT_HPP
#define GEOS_RESIDUALDUMPEVENT_HPP

#include "PeriodicEvent.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geos
{

class ResidualDumpEvent : public PeriodicEvent
{

public:

  ResidualDumpEvent( const string & name, Group * const parent );
  virtual ~ResidualDumpEvent() = default;
  static string catalogName() { return "ResidualDumpEvent"; };

  virtual void estimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    DomainPartition & domain ) override;


//           virtual void validate() const override;

  struct viewKeyStruct : PeriodicEvent::viewKeyStruct
  {
    static constexpr char const * secondaryTargetString() { return "secondaryTarget";}

    dataRepository::ViewKey secondaryTarget = { secondaryTargetString() };
  };

  void getTargetReferences() override;

  bool execute( real64 const time_n, real64 const dt, const integer cycleNumber,
                integer const,
                real64 const,
                DomainPartition & domain ) override;

  void validate() const override;


private:
  string m_secondaryEventTarget;
  ExecutableGroup * m_secondaryTarget;
  BitNodes< SolverBase::SolverGroupFlags > const * m_flagSecondaryTrigger;
};

} // geos

#endif //GEOS_RESIDUALDUMPEVENT_HPP
