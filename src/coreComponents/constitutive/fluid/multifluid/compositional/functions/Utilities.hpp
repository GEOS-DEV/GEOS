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
 * @file Utilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_UTILITIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_UTILITIES_HPP_

#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"

namespace geos
{
namespace constitutive
{
namespace compositional
{

/**
 * @brief Set a value to zero. Used for array initialisation
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static void setZero( real64 & value ){ value = 0.0; }

/**
 * @brief Normalise a composition in place to ensure that the components add up to unity
 * @param[in] numComps number of components
 * @param[in/out] composition composition to be normalized
 * @return the sum of the given values
 */
template< integer USD >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static real64 normalizeComposition( integer const numComps,
                                    arraySlice1d< real64, USD > const & composition )
{
  real64 totalMoles = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    totalMoles += composition[ic];
  }
  real64 const oneOverTotalMoles = 1.0 / (totalMoles + MultiFluidConstants::epsilon);
  for( integer ic = 0; ic < numComps; ++ic )
  {
    composition[ic] *= oneOverTotalMoles;
  }
  return totalMoles;
}

/**
 * @brief Checks if a composition array has numerically small values
 * @param[in] numComps number of values in the composition array
 * @param[in] composition composition to be checked
 * @return true if there is one zero value
 */
template< integer USD >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static bool hasZero( integer const numComps,
                     arraySlice1d< real64 const, USD > const & composition )
{
  for( integer ic = 0; ic < numComps; ++ic )
  {
    if( composition[ic] < MultiFluidConstants::epsilon )
    {
      return true;
    }
  }
  return false;
}

/**
 * @brief Calculate which components are present.
 * @details Creates a list of indices whose components have non-zero mole fraction.
 * @param[in] numComps number of components
 * @param[in] composition the composition of the fluid
 * @param[out] presentComponents the list of present components
 * @return the number of present components
 */
template< integer USD, typename ARRAY >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static integer calculatePresentComponents( integer const numComps,
                                           arraySlice1d< real64 const, USD > const & composition,
                                           ARRAY & presentComponents )
{
  // Check for machine-zero feed values
  integer presentCount = 0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    if( MultiFluidConstants::minForSpeciesPresence < composition[ic] )
    {
      presentComponents[presentCount++] = ic;
    }
  }
  presentComponents.resize( presentCount );
  return presentCount;
}

/**
 * @brief Determines whether a given composition represents a pure or nearly pure mixture.
 *
 * This function checks the mole fractions in a composition array to determine if the mixture
 * consists of a single dominant component (pure or almost pure). It identifies the index of the
 * component that is either strictly pure (all others are zero) or nearly pure (others are below a threshold).
 *
 * @tparam USD The stride or layout parameter for the array slice.
 *
 * @param[in] numComps Number of components in the composition array.
 * @param[in] composition Array slice containing the mole fractions of each component.
 * @param[out] pureComponent Index of the component that is strictly pure (only one non-zero mole fraction),
 *                           or -1 if no such component exists.
 * @param[out] almostPureComponent Index of the component that is nearly pure (others below a threshold),
 *                                 or -1 if no such component exists.
 * @return true if a strictly pure component is found (i.e., only one component has a significant mole fraction),
 *         false otherwise.
 */
template< integer USD >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static bool checkPureMixture( integer const numComps,
                              arraySlice1d< real64 const, USD > const & composition,
                              integer & pureComponent,
                              integer & almostPureComponent )
{
  pureComponent = -1;
  almostPureComponent = -1;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const zi1 = 1.0 - composition[ic];
    if( zi1 < MultiFluidConstants::minForSpeciesPresence )
    {
      pureComponent = ic;
    }
    if( zi1 < MultiFluidConstants::almostPureThreshold )
    {
      almostPureComponent = ic;
    }
  }
  return (0 <= pureComponent);
}

/**
 * @brief Solve the lineat system for the derivatives of the flash
 * @param[in/out] A the coefficient matrix. Destroyed after call
 * @param[in/out] X the rhs and solution
 * @return @c true if the problem is well solved @c false otherwise
 */
template< int USD >
GEOS_HOST_DEVICE
static bool solveLinearSystem( arraySlice2d< real64, USD > const & A,
                               arraySlice2d< real64, USD > const & X )
{
  // The flash and stability-test systems are A: N x N and X: N x (N+1).
  // denseLinearAlgebra::solve is GEOS_HOST_DEVICE but needs both dimensions at compile time.
  switch( A.size( 0 ) )
  {
    case 2:  return denseLinearAlgebra::solve< 2, 3 >( A, X );
    case 3:  return denseLinearAlgebra::solve< 3, 4 >( A, X );
    case 4:  return denseLinearAlgebra::solve< 4, 5 >( A, X );
    case 5:  return denseLinearAlgebra::solve< 5, 6 >( A, X );
    case 6:  return denseLinearAlgebra::solve< 6, 7 >( A, X );
    case 7:  return denseLinearAlgebra::solve< 7, 8 >( A, X );
    case 8:  return denseLinearAlgebra::solve< 8, 9 >( A, X );
    case 9:  return denseLinearAlgebra::solve< 9, 10 >( A, X );
    case 10:  return denseLinearAlgebra::solve< 10, 11 >( A, X );
    default: return false;
  }
}

} // namespace compositional
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_UTILITIES_HPP_
