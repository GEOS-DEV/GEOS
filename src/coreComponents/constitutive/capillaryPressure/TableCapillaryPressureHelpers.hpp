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
 * @file TableCapillaryPressureHelpers.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP

#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

struct TableCapillaryPressureHelpers
{

  /**
   * @brief Validate the capillary pressure table provided in input (increasing phase vol frac and cap pressure, etc)
   * @param[in] capPresTable the capillary pressure table (pc vs s) for a given phase)
   * @param[in] fullConstitutiveName the name of the constitutive model, for reporting if there is an issue
   * @param[in] capPresMustBeIncreasing flag saying that we expect an increasing cap pressure (otherwise, we expect a decreasing cap
   * pressure)
   */
  static
  void validateCapillaryPressureTable( TableFunction const & capPresTable,
                                       string const & fullConstitutiveName,
                                       bool const capPresMustBeIncreasing );

  static
  void validateCapillaryPressureTable( TableFunction const & capPresTable,
                                       string const & fullConstitutiveName,
                                       bool const capPresMustBeIncreasing,
                                       real64 & phaseMax,
                                       real64 & phaseMin );

  /**
   * @brief Populates the minimum phase volume fraction for each phase from the ends of the provided tables
   * @param[in] phaseOrder The ordering of the phases
   * @param[in] capPresTable The capillary pressure table
   * @param[out] minPhaseVolumeFraction The list of minimum phase volume fractions for each phase
   */
  static void populateMinPhaseVolumeFraction( arraySlice1d< integer const > const phaseOrder,
                                              TableFunction const & capPresTable,
                                              arraySlice1d< real64 > minPhaseVolumeFraction );

  /**
   * @brief Populates the minimum phase volume fraction for each phase from the ends of the provided tables
   * @param[in] phaseOrder The ordering of the phases
   * @param[in] capPresTableWettingIntermediate The wetting-intermediate capillary pressure table
   * @param[in] capPresTableNonWettingIntermediate The non-wetting-intermediate capillary pressure table
   * @param[out] minPhaseVolumeFraction The list of minimum phase volume fractions for each phase
   */
  static void populateMinPhaseVolumeFraction( arraySlice1d< integer const > const phaseOrder,
                                              TableFunction const & capPresTableWettingIntermediate,
                                              TableFunction const & capPresTableNonWettingIntermediate,
                                              arraySlice1d< real64 > minPhaseVolumeFraction );
};

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
