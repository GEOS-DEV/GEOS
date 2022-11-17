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
 * @file TableCapillaryPressureHelpers.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP

#include "functions/TableFunction.hpp"

namespace geosx
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
                                         real64& phaseMax,
                                         real64& phaseMin );

};

} // namespace constitutive

} // namespace geosx

#endif // GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
