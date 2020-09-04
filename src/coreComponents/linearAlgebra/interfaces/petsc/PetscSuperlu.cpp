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

namespace
{
/**
 * @brief Converts from GEOSX to SuperLU_Dist columns permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
string const & getColPermType( LinearSolverParameters::Direct::ColPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::ColPerm, string > const optionMap =
  {
    { LinearSolverParameters::Direct::ColPerm::none, "NATURAL" },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtplusA, "MMD_AT_PLUS_A" },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtA, "MMD_ATA" },
    { LinearSolverParameters::Direct::ColPerm::colAMD, "COLAMD" },
    { LinearSolverParameters::Direct::ColPerm::metis, "METIS_AT_PLUS_A" },
    { LinearSolverParameters::Direct::ColPerm::parmetis, "PARMETIS" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return optionMap.at( value );
}

/**
 * @brief Converts from GEOSX to SuperLU_Dist rows permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
string const & getRowPermType( LinearSolverParameters::Direct::RowPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::RowPerm, string > const optionMap =
  {
    { LinearSolverParameters::Direct::RowPerm::none, "NOROWPERM" },
    { LinearSolverParameters::Direct::RowPerm::mc64, "LargeDiag_MC64" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return optionMap.at( value );
}
}

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
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_colperm", getColPermType( params.direct.colPerm ).c_str() );
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_rowperm", getRowPermType( params.direct.rowPerm ).c_str() );
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

  // Create process grid: the target is to have the process grid as square as possible
  int const num_procs = MpiWrapper::Comm_size( matrix.getComm() );
  int prows = static_cast< int >( std::sqrt( num_procs ) );
  while( num_procs % prows )
  {
    --prows;
  }
  int pcols = num_procs/prows;
  std::tie( prows, pcols ) = std::minmax( prows, pcols );

  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_r", std::to_string( prows ).c_str() );
  PetscOptionsSetValue( nullptr, "-mat_superlu_dist_c", std::to_string( pcols ).c_str() );
}

}
