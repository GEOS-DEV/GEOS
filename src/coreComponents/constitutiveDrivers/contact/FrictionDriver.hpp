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

#ifndef GEOS_frictionDRIVER_HPP_
#define GEOS_frictionDRIVER_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{

class FrictionDriver : public TaskBase
{

public:
  FrictionDriver( const string & name,
                  Group * const parent );

  static string catalogName()
  { return "FrictionDriver"; }

  void postInputInitialization() override;

  virtual bool execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                        real64 const GEOS_UNUSED_PARAM( dt ),
                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition &
                        GEOS_UNUSED_PARAM( domain ) ) override;

  // /**
  //  * @brief Run test using loading protocol in table
  //  * @param i friction constitutive model
  //  * @param table Table with input/output time history
  //  */
  // template< typename friction_TYPE >
  // std::enable_if_t< std::is_same< constitutive::TableRelativePermeabilityHysteresis, friction_TYPE >::value, void >
  // runTest( friction_TYPE & friction,
  //          const arrayView2d< real64, 1 > & table );

  template< typename FRICTION_TYPE >
  void
  runTest( FRICTION_TYPE & friction,
           const arrayView2d< real64, 1 > & table );

  /**
   * @brief Ouput table to file for easy plotting
   */
  void outputResults();

  /**
   * @brief Read in a baseline table from file and compare with computed one (for unit testing purposes)
   */
  void compareWithBaseline();

private:

  template< typename FRICTION_TYPE >
  void resizeTable();

  template< typename FRICTION_TYPE >
  void resizeTables();

  // template< typename friction_TYPE >
  // std::enable_if_t< !std::is_same< constitutive::TableRelativePermeabilityHysteresis, friction_TYPE >::value, void >
  // resizeTable();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * frictionNameString()
    { return "friction"; }

    constexpr static char const * numStepsString()
    { return "steps"; }

    constexpr static char const * jumpFunctionString()
    { return "jumpControl"; }

    constexpr static char const * tractionFunctionString()
    { return "tractionControl"; }

    constexpr static char const * thetaString()
    { return "xTiltAngle";}

    constexpr static char const * phiString()
    { return "yTiltAngle";}

    constexpr static char const * outputString()
    { return "output"; }

    constexpr static char const * baselineString()
    { return "baseline"; }

  };

  integer m_numSteps;      ///< Number of load steps
  static integer const m_numColumns = 9;    ///< Number of columns in dat
  enum columnKeys { TIME, NJUMP, SLIP0, SLIP1, NTRAC, STRAC0, STRAC1, FS, TLIM };

  string m_jumpFunctionName; ///<
  string m_tractionFunctionName; ///<

  float m_theta, m_phi;///< x- and y-tilt of fault

  string m_frictionName;               ///< frictionType identifier
  string m_outputFile;              ///< Output file (optional, no output if not specified)

  array2d< real64 > m_table; ///< Table storing time-history of input/output

  Path m_baselineFile; ///< Baseline file (optional, for unit testing of solid models)


  static constexpr real64 m_baselineTol = 1e-3; ///< Comparison tolerance for baseline results
};


}

#endif //GEOS_FRICTIONDRIVER_HPP_
