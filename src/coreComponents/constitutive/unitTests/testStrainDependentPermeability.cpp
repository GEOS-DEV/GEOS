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

// Source includes
#include "constitutive/permeability/StrainDependentPermeability.hpp"

// TPL includes
#include <gtest/gtest.h>

#include <cmath>

using namespace geos;
using namespace geos::constitutive;

namespace
{

/// Initial porosity used in the tests
real64 constexpr phi0 = 0.2;

/// Minimum porosity of the model (see StrainDependentPermeabilityUpdate)
real64 constexpr minPorosity = 1.0e-6;

/// Analytical permeability (Chin et al., 2000)
real64 analyticalPermeability( real64 const k0,
                               real64 const initialPorosity,
                               real64 const n,
                               real64 const volStrain )
{
  real64 const phi = 1.0 - ( 1.0 - initialPorosity ) * std::exp( -volStrain );
  return k0 * std::pow( phi / initialPorosity, n );
}

class StrainDependentPermeabilityTest : public ::testing::Test
{
protected:

  static constexpr localIndex numElems = 1;
  static constexpr localIndex numQuad = 2;

  StrainDependentPermeabilityTest()
  {
    m_permeability.resize( numElems, 1, 3 );
    m_dPerm_dPressure.resize( numElems, 1, 3 );
    m_dPerm_dVolStrain.resize( numElems, 1, 3 );
    m_referencePermeability.resize( numElems, 1, 3 );
    m_volStrain.resize( numElems, numQuad );
    m_volStrain_n.resize( numElems, numQuad );

    for( integer dim = 0; dim < 3; ++dim )
    {
      m_referencePermeability[0][0][dim] = m_k0[dim];
      m_permeability[0][0][dim] = m_k0[dim];
      m_exponents[dim] = m_n[dim];
    }
  }

  StrainDependentPermeabilityUpdate createKernelWrapper()
  {
    return StrainDependentPermeabilityUpdate( m_permeability.toView(),
                                              m_dPerm_dPressure.toView(),
                                              m_dPerm_dVolStrain.toView(),
                                              m_referencePermeability.toViewConst(),
                                              m_volStrain.toView(),
                                              m_volStrain_n.toViewConst(),
                                              m_exponents );
  }

  /// Compute the permeability and its derivative for a given volumetric strain
  void compute( real64 const initialPorosity, real64 const volStrain )
  {
    StrainDependentPermeabilityUpdate const cmw = createKernelWrapper();
    cmw.compute( initialPorosity,
                 m_referencePermeability.toViewConst()[0][0],
                 volStrain,
                 m_permeability[0][0],
                 m_dPerm_dVolStrain[0][0] );
  }

  real64 const m_k0[3] = { 1.0e-15, 2.0e-15, 5.0e-16 };
  real64 const m_n[3] = { 3.0, 3.0, 2.0 };

  R1Tensor m_exponents;

  array3d< real64 > m_permeability;
  array3d< real64 > m_dPerm_dPressure;
  array3d< real64 > m_dPerm_dVolStrain;
  array3d< real64 > m_referencePermeability;
  array2d< real64 > m_volStrain;
  array2d< real64 > m_volStrain_n;
};

TEST_F( StrainDependentPermeabilityTest, zeroStrain )
{
  compute( phi0, 0.0 );

  for( integer dim = 0; dim < 3; ++dim )
  {
    EXPECT_NEAR( m_permeability[0][0][dim], m_k0[dim], 1.0e-12 * m_k0[dim] );
    // dk/de_v = n k (1 - phi_0) / phi_0 at zero strain
    real64 const expectedDerivative = m_n[dim] * m_k0[dim] * ( 1.0 - phi0 ) / phi0;
    EXPECT_NEAR( m_dPerm_dVolStrain[0][0][dim], expectedDerivative, 1.0e-12 * expectedDerivative );
  }
}

TEST_F( StrainDependentPermeabilityTest, analyticalValues )
{
  // negative: compaction, positive: dilation (extension positive)
  for( real64 const volStrain : { -0.05, -0.01, 0.01, 0.05 } )
  {
    compute( phi0, volStrain );

    for( integer dim = 0; dim < 3; ++dim )
    {
      real64 const expected = analyticalPermeability( m_k0[dim], phi0, m_n[dim], volStrain );
      EXPECT_NEAR( m_permeability[0][0][dim], expected, 1.0e-12 * expected );

      // compaction reduces the permeability, dilation increases it
      if( volStrain < 0.0 )
      {
        EXPECT_LT( m_permeability[0][0][dim], m_k0[dim] );
      }
      else
      {
        EXPECT_GT( m_permeability[0][0][dim], m_k0[dim] );
      }
    }
  }
}

TEST_F( StrainDependentPermeabilityTest, derivative )
{
  real64 const eps = 1.0e-7;

  for( real64 const volStrain : { -0.05, 0.0, 0.05 } )
  {
    real64 permPlus[3];
    real64 permMinus[3];
    real64 derivative[3];

    compute( phi0, volStrain + eps );
    for( integer dim = 0; dim < 3; ++dim )
    {
      permPlus[dim] = m_permeability[0][0][dim];
    }

    compute( phi0, volStrain - eps );
    for( integer dim = 0; dim < 3; ++dim )
    {
      permMinus[dim] = m_permeability[0][0][dim];
    }

    compute( phi0, volStrain );
    for( integer dim = 0; dim < 3; ++dim )
    {
      derivative[dim] = m_dPerm_dVolStrain[0][0][dim];
      real64 const finiteDifference = ( permPlus[dim] - permMinus[dim] ) / ( 2.0 * eps );
      EXPECT_NEAR( derivative[dim], finiteDifference, 1.0e-6 * LvArray::math::abs( derivative[dim] ) );
    }
  }
}

TEST_F( StrainDependentPermeabilityTest, clampedPorosity )
{
  // strong compaction: (1 - phi_0) exp( -e_v ) > 1, the porosity would become negative
  real64 const volStrain = std::log( 1.0 - phi0 ) - 0.1;
  compute( phi0, volStrain );

  for( integer dim = 0; dim < 3; ++dim )
  {
    real64 const expected = m_k0[dim] * std::pow( minPorosity / phi0, m_n[dim] );
    EXPECT_TRUE( std::isfinite( m_permeability[0][0][dim] ) );
    EXPECT_NEAR( m_permeability[0][0][dim], expected, 1.0e-12 * expected );
    EXPECT_DOUBLE_EQ( m_dPerm_dVolStrain[0][0][dim], 0.0 );
  }
}

TEST_F( StrainDependentPermeabilityTest, zeroInitialPorosity )
{
  compute( 0.0, 0.05 );

  for( integer dim = 0; dim < 3; ++dim )
  {
    EXPECT_DOUBLE_EQ( m_permeability[0][0][dim], m_k0[dim] );
    EXPECT_DOUBLE_EQ( m_dPerm_dVolStrain[0][0][dim], 0.0 );
  }
}

TEST_F( StrainDependentPermeabilityTest, updateFromVolumetricStrainIncrement )
{
  StrainDependentPermeabilityUpdate const cmw = createKernelWrapper();

  // converged volumetric strain at the two quadrature points
  m_volStrain_n[0][0] = 0.01;
  m_volStrain_n[0][1] = 0.03;

  // first quadrature point: the strain is stored but the permeability is not updated yet
  cmw.updateFromVolumetricStrainIncrement( 0, 0, phi0, 0.0 );
  EXPECT_DOUBLE_EQ( m_volStrain[0][0], 0.01 );
  for( integer dim = 0; dim < 3; ++dim )
  {
    EXPECT_DOUBLE_EQ( m_permeability[0][0][dim], m_k0[dim] );
  }

  // last quadrature point: the permeability is computed from the element-averaged strain
  cmw.updateFromVolumetricStrainIncrement( 0, 1, phi0, 0.02 );
  EXPECT_DOUBLE_EQ( m_volStrain[0][1], 0.05 );

  real64 const avgVolStrain = 0.5 * ( 0.01 + 0.05 );
  for( integer dim = 0; dim < 3; ++dim )
  {
    real64 const expected = analyticalPermeability( m_k0[dim], phi0, m_n[dim], avgVolStrain );
    EXPECT_NEAR( m_permeability[0][0][dim], expected, 1.0e-12 * expected );
  }
}

} // namespace
