/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MixedVEMCompliance.cpp
 */

#include "mixedVEM/MixedVEMCompliance.hpp"

#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

#include <cmath>

namespace geos
{

namespace mixedVEM
{

void invertStiffness( real64 const (&stiffness)[NUM_SYM_COMP][NUM_SYM_COMP],
                      real64 (& compliance)[NUM_SYM_COMP][NUM_SYM_COMP] )
{
  array2d< real64 > source( NUM_SYM_COMP, NUM_SYM_COMP );
  array2d< real64 > target( NUM_SYM_COMP, NUM_SYM_COMP );

  for( integer i = 0; i < NUM_SYM_COMP; ++i )
  {
    for( integer j = 0; j < NUM_SYM_COMP; ++j )
    {
      source( i, j ) = stiffness[i][j];
    }
  }

  arrayView2d< real64 const > const sourceView = source.toViewConst();
  arrayView2d< real64 > const targetView = target.toView();

  BlasLapackLA::matrixInverse( sourceView.toSliceConst(), targetView.toSlice() );

  for( integer i = 0; i < NUM_SYM_COMP; ++i )
  {
    for( integer j = 0; j < NUM_SYM_COMP; ++j )
    {
      compliance[i][j] = target( i, j );
    }
  }
}

void convertVoigtStiffness( real64 const (&voigtStiffness)[NUM_SYM_COMP][NUM_SYM_COMP],
                            real64 (& stiffness)[NUM_SYM_COMP][NUM_SYM_COMP] )
{
  // sigma_o = S sigma_v and eps_v = S eps_o with S = diag(1,1,1,sqrt2,sqrt2,sqrt2),
  // so C_o = S C_v S
  real64 const s[NUM_SYM_COMP] = { 1.0, 1.0, 1.0, M_SQRT2, M_SQRT2, M_SQRT2 };

  for( integer i = 0; i < NUM_SYM_COMP; ++i )
  {
    for( integer j = 0; j < NUM_SYM_COMP; ++j )
    {
      stiffness[i][j] = s[i] * s[j] * voigtStiffness[i][j];
    }
  }
}

} // namespace mixedVEM

} // namespace geos
