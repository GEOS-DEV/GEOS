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
 * @file MixedVEMCellOutput.hpp
 *
 * Projection of the mixed VEM unknowns onto writable cell arrays.
 *
 * sigma_h is virtual inside E and its degrees of freedom sit on the faces, while u_h|_E
 * is a rigid body motion rather than a single vector, so neither unknown is a cell array
 * as it stands. The two computable cell representatives are Pi_E sigma_h, the constant
 * stress that the consistency term already reproduces exactly, and the value of u_h at
 * x_E, which is the translation part of the rigid motion. Both are what a finite element
 * run reports per cell, so they are the quantities to compare against.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMCELLOUTPUT_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMCELLOUTPUT_HPP_

#include "mixedVEM/MixedVEMElementOperators.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @brief Split the rigid motion coefficients into a displacement and a rotation.
 * @param[in] rigidMotion the six coefficients of u_h|_E on the basis {R_i}
 * @param[out] displacement u_h(x_E), the translation part
 * @param[out] rotation the infinitesimal rotation omega
 *
 * u_h(x) = alpha + omega ^ (x - x_E), so the value at the cell center is alpha.
 */
GEOS_HOST_DEVICE
inline void computeCellDisplacement( real64 const * const rigidMotion,
                                     real64 (& displacement)[3],
                                     real64 (& rotation)[3] )
{
  for( integer i = 0; i < 3; ++i )
  {
    displacement[i] = rigidMotion[i];
    rotation[i] = rigidMotion[3 + i];
  }
}

/**
 * @brief Evaluate u_h|_E at an arbitrary point of the cell.
 * @param[in] rigidMotion the six coefficients of u_h|_E on the basis {R_i}
 * @param[in] elemCenter the point x_E
 * @param[in] point the evaluation point
 * @param[out] value the displacement at @p point
 */
GEOS_HOST_DEVICE
inline void evaluateCellDisplacement( real64 const * const rigidMotion,
                                      real64 const (&elemCenter)[3],
                                      real64 const (&point)[3],
                                      real64 (& value)[3] )
{
  real64 r[3];
  for( integer i = 0; i < 3; ++i )
  {
    r[i] = point[i] - elemCenter[i];
  }

  real64 omega[3] = { rigidMotion[3], rigidMotion[4], rigidMotion[5] };
  LvArray::tensorOps::crossProduct( value, omega, r );

  for( integer i = 0; i < 3; ++i )
  {
    value[i] += rigidMotion[i];
  }
}

/**
 * @brief Project the stress degrees of freedom onto the constant cell stress.
 * @param[in] projection the matrix P_E
 * @param[in] numFaces the number of faces of the element
 * @param[in] stressDof the 6 numFaces stress degrees of freedom of the element
 * @param[out] stress Pi_E sigma_h in the order (xx, yy, zz, yz, xz, xy)
 *
 * P_E returns the coefficients on the orthonormal basis, where the shear entries carry a
 * factor sqrt(2); the output drops it so the six components are the plain tensor entries
 * and match the convention of the other solvers.
 */
GEOS_HOST_DEVICE
inline void computeCellStress( MatrixSliceConst const & projection,
                               integer const numFaces,
                               real64 const * const stressDof,
                               real64 (& stress)[NUM_SYM_COMP] )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    real64 sum = 0.0;
    for( integer j = 0; j < numStressDof; ++j )
    {
      sum += projection( a, j ) * stressDof[j];
    }
    stress[a] = ( a < 3 ) ? sum : INV_SQRT_2 * sum;
  }
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMCELLOUTPUT_HPP_
