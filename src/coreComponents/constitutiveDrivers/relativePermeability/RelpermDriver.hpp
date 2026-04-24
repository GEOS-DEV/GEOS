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

#ifndef GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVER_HPP
#define GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVER_HPP

#include "constitutiveDrivers/ConstitutiveDriver.hpp"

namespace geos
{

namespace constitutive
{
class RelativePermeabilityBase;
}

/**
 * @class RelpermDriver
 *
 * Class to allow for testing Relative permeability models without the
 * complexity of setting up a full simulation.
 */
class RelpermDriver : public ConstitutiveDriver
{
public:
  RelpermDriver( const string & name, Group * const parent );

  static string catalogName() { return "RelpermDriver"; }

  void postInputInitialization() override;

  bool execute() override;

  void getColumnNames( string_array & columnNames ) const override;

  /**
   * @brief Run test using loading protocol in table
   * @param relperm Relperm constitutive model
   * @param table Table with input/output time history
   */
  template< typename RELPERM_TYPE >
  void runTest( RELPERM_TYPE & relperm, const arrayView2d< real64, 1 > & table );

private:
  /**
   * @brief Get the relative permeability model from the catalog
   */
  constitutive::RelativePermeabilityBase & getRelperm();
  constitutive::RelativePermeabilityBase const & getRelperm() const;

  void allocateTable( integer numColumns, integer numPhases );

  /**
   * @brief Initialises the table by filling in primary variables
   */
  void initializeTable();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * relpermNameString() { return "relperm"; }
    constexpr static char const * historicalSaturationsString() { return "historicalSaturations"; }
  };

  string m_relpermName;               ///< relPermType identifier
  array1d< real64 > m_historicalSaturations;
};

}

#endif //GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVER_HPP
