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
 * @file EquilibriumReaction.cpp
 */

#include "constitutive/fluid/PVTFunctions/EquilibriumReaction.hpp"

#include "constitutive/fluid/PVTFunctions/CO2EOSSolver.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{


} // namespace

EquilibriumReaction::EquilibriumReaction( string const & name,
                                          string_array const & inputParams,
                                          string_array const & phaseNames,
                                          string_array const & componentNames,
                                          array1d< real64 > const & componentMolarWeight ):
  ReactionBase( name,
                componentNames,
                componentMolarWeight )
{
  GEOSX_THROW_IF_NE_MSG( phaseNames.size(), 2,
                         "The EquilibriumReaction model is a two-phase model",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( componentNames.size(), 2,
                         "The EquilibriumReaction model is a two-component model",
                         InputError );
}

EquilibriumReaction::KernelWrapper EquilibriumReaction::createKernelWrapper() const
{
  return KernelWrapper(  );
}

REGISTER_CATALOG_ENTRY( ReactionBase, EquilibriumReaction, string const &, string_array const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
