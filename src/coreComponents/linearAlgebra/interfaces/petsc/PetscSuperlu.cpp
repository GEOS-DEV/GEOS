/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscSuperlu.cpp
 */

#include "PetscSuperlu.hpp"
#include <petsc.h>
#include <superlu_ddefs.h>

namespace geosx
{

void SuperLU_DistSetFromOptions( PetscMatrix const & matrix,
                                 LinearSolverParameters const & params )
{
  // Initialize options.
  if( params.logLevel > 0 )
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_statprint", "1" );
  }
  else
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_statprint", "0" );
  }

  if( params.direct.equilibrate )
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_equil", "1" );
  }
  else
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_equil", "0" );
  }
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_colperm", PetscGetColPermType( params.direct.colPerm ).c_str() );
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_rowperm", PetscGetRowPermType( params.direct.rowPerm ).c_str() );
  if( params.direct.replaceTinyPivot )
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_replacetinypivot", "1" );
  }
  else
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_replacetinypivot", "0" );
  }
  if( params.direct.iterativeRefine )
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_iterrefine", "1" );
  }
  else
  {
    PetscOptionsSetValue( nullptr, "-mat_superlu_dist_iterrefine", "0" );
  }

  // Create process grid.
  int const num_procs = MpiWrapper::Comm_size( matrix.getComm() );
  int pcols = 1;
  int prows = 1;
  while( prows*pcols <= num_procs )
  {
    ++prows;
  }
  --prows;
  pcols = num_procs/prows;
  while( prows*pcols != num_procs )
  {
    prows -= 1;
    pcols = num_procs/prows;
  }
  char prowsChar[100] = { 0 };
  sprintf( prowsChar, "%i", prows );
  char pcolsChar[100] = { 0 };
  sprintf( pcolsChar, "%i", pcols );

  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_r", prowsChar );
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_c", pcolsChar );
}

string const & PetscGetColPermType( string const & value )
{
  static std::map< string, string > const optionMap =
  {
    { "none", "NATURAL" },
    { "MMD_At+A", "MMD_AT_PLUS_A" },
    { "MMD_AtA", "MMD_ATA" },
    { "metis", "METIS_AT_PLUS_A" },
    { "parmetis", "PARMETIS" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return optionMap.at( value );
}

string const & PetscGetRowPermType( string const & value )
{
  static std::map< string, string > const optionMap =
  {
    { "none", "NOROWPERM" },
    { "mc64", "LargeDiag_MC64" },
    { "awpm", "LargeDiag_AWPM" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return optionMap.at( value );
}

}
