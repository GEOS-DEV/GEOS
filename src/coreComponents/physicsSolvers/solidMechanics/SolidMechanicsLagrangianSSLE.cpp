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
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#include "SolidMechanicsLagrangianSSLE.hpp"

#include "codingUtilities/Utilities.hpp"
#include "finiteElement/Kinematics.h"
#include "discretizationMethods/NumericalMethodsManager.hpp"



namespace geosx
{

using namespace constitutive;

SolidMechanicsLagrangianSSLE::SolidMechanicsLagrangianSSLE( string const & name,
                                                            Group * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->m_strainTheory = 0;
}

SolidMechanicsLagrangianSSLE::~SolidMechanicsLagrangianSSLE()
{}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianSSLE, string const &, dataRepository::Group * const )
} /* namespace geosx */
