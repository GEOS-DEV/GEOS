/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscInterface.cpp
 */

#include "PetscInterface.hpp"

#include <petscsys.h>

namespace geosx
{

void PetscInterface::initialize( int & argc, char **& argv )
{
  // PetscOptionsSetValue( nullptr, "-log_view", "" );
  // PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr );
  if( argc > 0 )
  {
    PetscInitialize( &argc, &argv, nullptr, nullptr );
  }
  else
  {
    PetscInitializeNoArguments();
  }
  PetscPopSignalHandler();
}

void PetscInterface::finalize()
{
  PetscFinalize();
}

} //namespace geosx
