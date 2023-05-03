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
 * @file SolidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_

namespace geos
{
namespace constitutive
{

struct SolidUtilities
{

  /**
   * @brief Perform a finite-difference check of the stiffness computation
   *
   * This method uses several stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix, and then computes an error between
   * the coded stiffness method and the finite difference version.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness.
   *
   * @param solid the solid kernel wrapper
   * @param k the element number
   * @param q the quadrature index
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   * @param print flag to decide if debug output is printed or not
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static bool
  checkSmallStrainStiffness( SOLID_TYPE const & solid,
                             localIndex k,
                             localIndex q,
                             real64 const ( &strainIncrement )[6],
                             bool print = false )
  {
    real64 stiffness[6][6]{};     // coded stiffness
    real64 stiffnessFD[6][6]{};   // finite difference approximation
    real64 stress[6]{};           // original stress

    solid.smallStrainUpdate( k, q, strainIncrement, stress, stiffness );
    SolidUtilities::computeSmallStrainFiniteDifferenceStiffness( solid, k, q, strainIncrement, stiffnessFD );

    // compute relative error between two versions

    real64 error = 0;
    real64 norm = 0;

    for( localIndex i=0; i<6; ++i )
    {
      for( localIndex j=0; j<6; ++j )
      {
        error += fabs( stiffnessFD[i][j]-stiffness[i][j] );
        norm += fabs( stiffnessFD[i][j] );
      }
    }
    error /= norm;

    // optional printing for debugging purposes

    if( print )
    {
      for( localIndex i=0; i<6; ++i )
      {
        for( localIndex j=0; j<6; ++j )
        {
          printf( "[%8.1e vs %8.1e] ", stiffnessFD[i][j], stiffness[i][j] );
        }
        printf( "\n" );
      }
    }

    return (error < 1e-3);
  }

  /**
   * @brief Perform a finite-difference stiffness computation
   *
   * This method uses stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness, as it will primarily be used to check
   * the hand coded tangent against a finite difference reference.
   * A similar method would need to be implemented to check compressed stiffness,
   * stress-only, or finite-strain interfaces.
   *
   * @param solid the solid kernel wrapper
   * @param k the element number
   * @param q the quadrature index
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   * @param stiffnessFD finite different stiffness approximation
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  computeSmallStrainFiniteDifferenceStiffness( SOLID_TYPE const & solid,
                                               localIndex k,
                                               localIndex q,
                                               real64 const ( &strainIncrement )[6],
                                               real64 ( & stiffnessFD )[6][6] )
  {
    real64 stiffness[6][6]{};      // coded stiffness
    real64 stress[6]{};            // original stress
    real64 stressFD[6]{};          // perturbed stress
    real64 strainIncrementFD[6]{}; // perturbed strain
    real64 norm = 0;             // norm for scaling (note: method is fragile w.r.t. scaling)

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] = strainIncrement[i];
      norm += fabs( strainIncrement[i] );
    }

    real64 eps = 1e-4*norm; // finite difference perturbation

    solid.smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] += eps;

      if( i>0 )
      {
        strainIncrementFD[i-1] -= eps;
      }

      solid.smallStrainUpdate( k, q, strainIncrementFD, stressFD, stiffnessFD );

      for( localIndex j=0; j<6; ++j )
      {
        stiffnessFD[j][i] = (stressFD[j]-stress[j])/eps;
      }
    }
  }

  /**
   * @brief Hypo update (small strain, large rotation).
   *
   * This function uses a call to the small strain update, followed by
   * a rotation correction using the Hughes-Winget incrementally objective
   * algorithm.  One can imagine the material deforming in small strain, but in
   * a reference frame that rotates to track the body's local rotation.  This
   * provides a convenient way to extend small strain models to a finite rotation
   * regime.  From the assumption of small deformations, the Cauchy stress and
   * Kirchoff stress are approximately equal (det F ~ 1).
   *
   * We use a post-rotation of the new stress (as opposed to pre-rotation of the old stress)
   * as we don't have to unrotate the old stress in the event a rewind is required.
   * We do not post-rotate the stiffness tensor, which is an approximation, but
   * should be sufficient for small rotation increments.
   *
   * This method does not work for anisotropic properties or yield functions
   * (without some care) and a co-rotational formulation should be considered instead.
   *
   * @param solid the solid kernel wrapper
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
   * @param[in] Rot The incremental rotation tensor
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate( SOLID_TYPE const & solid,
              localIndex const k,
              localIndex const q,
              real64 const ( &Ddt )[6],
              real64 const ( &Rot )[3][3],
              real64 ( & stress )[6],
              real64 ( & stiffness )[6][6] )
  {
    solid.smallStrainUpdate( k, q, Ddt, stress, stiffness );

    real64 temp[6]{};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, stress );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }

/**
 * @brief Hypo update, returning only stress
 *
 * @param solid the solid kernel wrapper
 * @param[in] k The element index.
 * @param[in] q The quadrature point index.
 * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
 * @param[in] Rot The incremental rotation tensor
 * @param[out] stress New stress value (Cauchy stress)
 */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate_StressOnly( SOLID_TYPE const & solid,
                         localIndex const k,
                         localIndex const q,
                         real64 const ( &Ddt )[6],
                         real64 const ( &Rot )[3][3],
                         real64 ( & stress )[6] )
  {
    solid.smallStrainUpdate_StressOnly( k, q, Ddt, stress );

    real64 temp[6]{};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, stress );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }


};

} // namespace constitutive

} // namespace geos

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_ */
