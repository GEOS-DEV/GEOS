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


// Source includes
// #include "mainInterface/initialization.hpp"
// #include "finiteVolume/CellElementStencilTPFA.hpp"
// #include "physicsSolvers/fluidFlow/kernels/singlePhase/unitTests/testFlowKernelHelpers.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

// TPL includes
#include <cmath>
#include <gtest/gtest.h>


using namespace geos;
constexpr real64 EPS = 1e-16;


TEST(testIsoTAligned, Velocity_Tetra)
{

  constexpr int NFACES = 4;    //2-hexa problem 
  constexpr int NTESTS = 2;
  real64 tests[NTESTS][NFACES] = {{1.,1.,1.,1.},
 {2.,-1.,2.,-1.}};
 real64 expectedVel[NTESTS][3] = {
  {0.,0.,0.},
  {0.,3.*LvArray::math::sqrt(3)/2.,0.}
 };
  real64 const normals[NFACES][3] =  {
    {1./LvArray::math::sqrt(3), 1./LvArray::math::sqrt(3), 1./LvArray::math::sqrt(3)},
    {1./LvArray::math::sqrt(3), -1./LvArray::math::sqrt(3), -1./LvArray::math::sqrt(3)},
    {-1./LvArray::math::sqrt(3), 1./LvArray::math::sqrt(3), -1./LvArray::math::sqrt(3)},
    {-1./LvArray::math::sqrt(3), -1./LvArray::math::sqrt(3), 1./LvArray::math::sqrt(3)}
  };


for(int itest=0; itest<NTESTS; itest++){
  array1d< real64 > fluxes; 
  fluxes.resize(NFACES);
  fluxes[0] = tests[itest][0];
  fluxes[1] = tests[itest][1];
  fluxes[2] = tests[itest][2];
  fluxes[3] = tests[itest][3];

  auto phaseVelocity = computationalGeometry::computeVelocities_<NFACES>(normals, fluxes.toViewConst());
  GEOS_LOG(GEOS_FMT("phaseVelocity {}", phaseVelocity));
  EXPECT_LT(phaseVelocity[0][0]-expectedVel[itest][0], EPS);
  EXPECT_LT(phaseVelocity[0][1]-expectedVel[itest][1], EPS);
  EXPECT_LT(phaseVelocity[0][2]-expectedVel[itest][2], EPS);

  }


}

TEST( testRotated3D, Velocity_Hexa )
{
  constexpr int NFACES = 6;    //2-hexa problem 
  constexpr int NTESTS = 2;
  real64 tests[NTESTS][NFACES] = {{1.,1.,1.,1.,1.,1.},
 {2.,-1.,2.,-1.,2.,-1.}};
 real64 expectedVel[NTESTS][3] = {
  {0.,0.,0.},
  {0.,1.5*LvArray::math::sqrt(2),1.5}
 };
  real64 const normals[NFACES][3] =  {{1./LvArray::math::sqrt(2), 1./LvArray::math::sqrt(2), 0.},
    {-1./LvArray::math::sqrt(2), -1./LvArray::math::sqrt(2), 0.},
    {-1./LvArray::math::sqrt(2), 1./LvArray::math::sqrt(2), 0.},
    {1./LvArray::math::sqrt(2), -1./LvArray::math::sqrt(2), 0.},
    {0., 0., 1.},
    {0., 0., -1.}
  };

for(int itest=0; itest<NTESTS; itest++){
  array1d< real64 > fluxes; 
  fluxes.resize(NFACES);
  //{1., 1., 1., 1., 1., 1.};
  fluxes[0] = tests[itest][0];
  fluxes[1] = tests[itest][1];
  fluxes[2] = tests[itest][2];
  fluxes[3] = tests[itest][3];
  fluxes[4] = tests[itest][4];
  fluxes[5] = tests[itest][5];

  // real64 const fluxes[nfaces] = 

  auto phaseVelocity = computationalGeometry::computeVelocities_<NFACES>(normals, fluxes.toViewConst());
  GEOS_LOG(GEOS_FMT("phaseVelocity {}", phaseVelocity));
  EXPECT_LT(phaseVelocity[0][0]-expectedVel[itest][0], EPS);
  EXPECT_LT(phaseVelocity[0][1]-expectedVel[itest][1], EPS);
  EXPECT_LT(phaseVelocity[0][2]-expectedVel[itest][2], EPS);

  }

}

TEST( testAligned3D, Velocity_Hexa )
{
  constexpr int NFACES = 6;    //2-hexa problem 
  constexpr int NTESTS = 2;
  real64 tests[NTESTS][NFACES] = {{1.,1.,1.,1.,1.,1.},
 {2.,-1.,2.,-1.,2.,-1.}};
 real64 expectedVel[NTESTS][3] = {
  {0.,0.,0.},
  {1.5,1.5,1.5}
 };

  real64 const normals[NFACES][3] =
  {{1., 0., 0.},
    {-1., 0., 0.},
    {0., 1., 0.},
    {0., -1., 0.},
    {0., 0., 1.},
    {0., 0., -1.}
  };//ROT 0


  for(int itest=0; itest<NTESTS; itest++){
  array1d< real64 > fluxes; 
  fluxes.resize(NFACES);
  //{1., 1., 1., 1., 1., 1.};
  fluxes[0] = tests[itest][0];
  fluxes[1] = tests[itest][1];
  fluxes[2] = tests[itest][2];
  fluxes[3] = tests[itest][3];
  fluxes[4] = tests[itest][4];
  fluxes[5] = tests[itest][5];

  // real64 const fluxes[nfaces] = 

  auto phaseVelocity = computationalGeometry::computeVelocities_<NFACES>(normals, fluxes.toViewConst());
  GEOS_LOG(GEOS_FMT("phaseVelocity {}", phaseVelocity));
  EXPECT_LT(phaseVelocity[0][0]-expectedVel[itest][0], EPS);
  EXPECT_LT(phaseVelocity[0][1]-expectedVel[itest][1], EPS);
  EXPECT_LT(phaseVelocity[0][2]-expectedVel[itest][2], EPS);

  }

}
