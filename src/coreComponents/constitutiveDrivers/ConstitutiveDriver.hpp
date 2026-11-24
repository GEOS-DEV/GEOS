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
 * @file ConstitutiveDriver.hpp
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_CONSTITUTIVEDRIVER_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_CONSTITUTIVEDRIVER_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{
namespace constitutive
{
class ConstitutiveManager;
}

/**
 * @class ConstitutiveDriver
 *
 * Base class for a constitutive driver.
 *
 */
class ConstitutiveDriver : public TaskBase
{
public:
  ConstitutiveDriver( const string & name, Group * const parent );
  ~ConstitutiveDriver() override = default;

  void postInputInitialization() override;

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  virtual bool execute() = 0;

protected:
  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * numStepsString() { return "steps"; }
    constexpr static char const * outputString() { return "output"; }
    constexpr static char const * precisionString() { return "precision"; }
    constexpr static char const * outputDataString() { return "outputData"; }
  };

  /**
   * @brief Validate results by checking residual and removing erroneous data
   */
  virtual bool validateResults();

  /**
   * @brief Ouput table to file or screen
   */
  void outputResults() const;

  /**
   * Get the column names
   */
  virtual void getColumnNames( string_array & columnNames ) const = 0;

  /**
   * Allocate space for the table
   */
  void allocateTable( integer const numColumns, real64 const minTime, real64 const maxTime );

  /**
   * Extract the constitutive manager
   */
  constitutive::ConstitutiveManager & getConstitutiveManager();
  constitutive::ConstitutiveManager const & getConstitutiveManager() const;

private:
  /**
   * Output the table to a file
   */
  void outputToFile() const;

  /**
   * Output the table to the console
   */
  void outputToConsole() const;

protected:
  static constexpr integer minPrecision = 1;
  static constexpr integer maxPrecision = 8;

  integer m_numSteps;           ///< Number of steps for experiment
  string m_outputFile;          ///< Output file (optional, output to screen if not specified)
  integer m_precision{4};       ///< The precision to use to data out to files

  array2d< real64 > m_table;    ///< Table storing time-history of data

  enum columnKeys { TIME };     ///< Enumeration of "input" column keys for readability
};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVEDRIVERS_CONSTITUTIVEDRIVER_HPP_ */
