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
 * @file ReactiveFluidDriver.hpp
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_FLUID_PVTDRIVER_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_FLUID_PVTDRIVER_HPP_

#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "events/tasks/TaskBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

using namespace constitutive;

/**
 * @class ReactiveFluidDriver
 *
 * Class to allow for testing PVT behavior without the
 * complexity of setting up a full simulation.
 *
 */
class ReactiveFluidDriver : public TaskBase
{
public:
  ReactiveFluidDriver( const string & name,
             Group * const parent );
  ~ReactiveFluidDriver() override;

  static string catalogName() { return "ReactiveFluidDriver"; }

  void postProcessInput() override;

  virtual bool execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                        real64 const GEOSX_UNUSED_PARAM( dt ),
                        integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override;

  /**
   * @brief Run test using loading protocol in table
   * @param fluid Fluid constitutive model
   * @param table Table with input/output time history
   */
  template< typename FLUID_TYPE >
  void runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table );

  /**
   * @brief Ouput table to file for easy plotting
   */
  void outputResults();

  /**
   * @brief Read in a baseline table from file and compare with computed one (for unit testing purposes)
   */
  void compareWithBaseline();

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * fluidNameString() { return "fluid"; }
    constexpr static char const * pressureFunctionString() { return "pressureControl"; }
    constexpr static char const * temperatureFunctionString() { return "temperatureControl"; }
    constexpr static char const * numStepsString() { return "steps"; }
    constexpr static char const * outputString() { return "output"; }
    constexpr static char const * baselineString() { return "baseline"; }
    constexpr static char const * feedString() { return "feedComposition"; }
  };

  integer m_numSteps;      ///< Number of load steps
  integer m_numColumns;    ///< Number of columns in data table (depends on number of fluid phases)
  integer m_numPhases;     ///< Number of fluid phases

  integer m_numPrimarySpecies; ///< Number of primary species (components)
  integer m_numSecondarySpecies; ///< Number of secondary species (components)

  string m_fluidName;               ///< Fluid identifier
  string m_pressureFunctionName;    ///< Time-dependent function controlling pressure
  string m_temperatureFunctionName; ///< Time-dependent function controlling temperature
  string m_outputFile;              ///< Output file (optional, no output if not specified)

  array1d< real64 > m_feed;  ///< User specified feed composition
  array2d< real64 > m_table; ///< Table storing time-history of input/output

  Path m_baselineFile; ///< Baseline file (optional, for unit testing of solid models)

  enum columnKeys { TIME, PRES, TEMP }; ///< Enumeration of "input" column keys for readability

  static constexpr real64 m_baselineTol = 1e-3; ///< Comparison tolerance for baseline results
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_FLUID_PVTDRIVER_HPP_ */
