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
 * @file TableRelativePermeabilityHelpers.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_TABLERELATIVEPERMEABILITYHELPERS_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_TABLERELATIVEPERMEABILITYHELPERS_HPP

#include "functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

struct TableRelativePermeabilityHelpers
{

  /**
   * @brief Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTable the relative permeability table (kr vs s) for a given phase)
   * @param[in] fullConstitutiveName name of the constitutive model
   * @param[out] phaseMinVolFrac the phase minimum volume fraction read from the table
   * @param[out] phaseMaxVolFrac the phase maximum volume fraction read from the table
   * @param[out] phaseRelPermMinEndPoint the end-point relative permeability
   */
  static
  void validateRelativePermeabilityTable( TableFunction const & relPermTable, string const & fullConstitutiveName,
                                          real64 & phaseMinVolFrac, real64 & phaseMaxVolFrac,
                                          real64 & phaseRelPermMinEndPoint, real64 & phaseRelPermMaxEndPoint );

};

} // namespace constitutive

} // namespace geosx

#endif // GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
