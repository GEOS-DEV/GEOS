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
 * @file CO2SolubilitySpycherPruess.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYSPYCHERPRUESS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYSPYCHERPRUESS_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

class TableFunction;
class FunctionManager;

namespace constitutive
{
namespace PVTProps
{

struct CO2SolubilitySpycherPruess
{

/**
 * @brief Create CO2 and H2O solubility table based on Spycher, Pruess, Ennis-King (2003)
 * @details The generated table is a 2D table with lookup properties pressure (in Pa) and
 *          temperature (in degC). The returned CO2 solubility is in mole of CO2 per kg of
 *          H2O and the returned water vapourisation is in moles of H2O per kg of CO2.
 * @param[in] inputParams A list of input parameters
 * @param[in] functionName The name of the model
 * @param[in] functionManager The function manager to which the table should be attached
 * @return The created tables with first CO2 solubility table and second H2O solubility table
 */
  static std::pair< TableFunction const *, TableFunction const * >
  makeSolubilityTables( string_array const & inputParams,
                        string const & functionName,
                        FunctionManager & functionManager );

};

} // end namespace PVTProps
} // end namespace constitutive
} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYSPYCHERPRUESS_HPP_
