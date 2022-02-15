/*	
 * ------------------------------------------------------------------------------------------------------------	
 * SPDX-License-Identifier: LGPL-2.1-only	
 *	
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC	
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University	
 * Copyright (c) 2018-2020 Total, S.A	
 * Copyright (c) 2020-     GEOSX Contributors	
 * All right reserved	
 *	
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.	
 * ------------------------------------------------------------------------------------------------------------	
 */

#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#include <gtest/gtest.h>

using namespace geosx;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( LinearSolverParametersEnums, SolverType )
{
  using EnumType = LinearSolverParameters::SolverType;

  ASSERT_EQ( "direct", toString( EnumType::direct ) );
  ASSERT_EQ( "cg", toString( EnumType::cg ) );
  ASSERT_EQ( "gmres", toString( EnumType::gmres ) );
  ASSERT_EQ( "fgmres", toString( EnumType::fgmres ) );
  ASSERT_EQ( "bicgstab", toString( EnumType::bicgstab ) );
  ASSERT_EQ( "preconditioner", toString( EnumType::preconditioner ) );
}


TEST( LinearSolverParametersEnums, PreconditionerType )
{
  using EnumType = LinearSolverParameters::PreconditionerType;

  ASSERT_EQ( "none", toString( EnumType::none ) );
  ASSERT_EQ( "jacobi", toString( EnumType::jacobi ) );
  ASSERT_EQ( "l1jacobi", toString( EnumType::l1jacobi ) );
  ASSERT_EQ( "fgs", toString( EnumType::fgs ) );
  ASSERT_EQ( "sgs", toString( EnumType::sgs ) );
  ASSERT_EQ( "l1sgs", toString( EnumType::l1sgs ) );
  ASSERT_EQ( "chebyshev", toString( EnumType::chebyshev ) );
  ASSERT_EQ( "iluk", toString( EnumType::iluk ) );
  ASSERT_EQ( "ilut", toString( EnumType::ilut ) );
  ASSERT_EQ( "ic", toString( EnumType::ic ) );
  ASSERT_EQ( "ict", toString( EnumType::ict ) );
  ASSERT_EQ( "amg", toString( EnumType::amg ) );
  ASSERT_EQ( "mgr", toString( EnumType::mgr ) );
  ASSERT_EQ( "block", toString( EnumType::block ) );
  ASSERT_EQ( "direct", toString( EnumType::direct ) );
  ASSERT_EQ( "bgs", toString( EnumType::bgs ) );
}


int main( int argc, char * * argv )
{
  geosx::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
