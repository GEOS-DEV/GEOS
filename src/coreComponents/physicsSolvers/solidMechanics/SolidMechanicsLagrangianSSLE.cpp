/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#include "SolidMechanicsLagrangianSSLE.hpp"

#include "codingUtilities/Utilities.hpp"
#include "finiteElement/Kinematics.h"



namespace geosx
{

using namespace constitutive;

SolidMechanicsLagrangianSSLE::SolidMechanicsLagrangianSSLE( string const & name,
                                                            ManagedGroup * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->m_strainTheory = 0;
}

SolidMechanicsLagrangianSSLE::~SolidMechanicsLagrangianSSLE()
{}




REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianSSLE, string const &, dataRepository::ManagedGroup * const )
} /* namespace geosx */

