/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC*
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FRICTIONDRIVERRUNTEST_HPP_
#define GEOS_FRICTIONDRIVERRUNTEST_HPP_

#include "constitutiveDrivers/contact/FrictionDriver.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"
#include "constitutive/solid/SolidFields.hpp"

#include <cmath>
#include <conduit.hpp>

namespace geos
{

template< typename FRICTION_TYPE >
void
FrictionDriver::runTest( FRICTION_TYPE & friction,
                         const arrayView2d< real64 > & table )
{

  array1d< integer > const ghostRank( 1 ); ghostRank[0] = -1;

  array1d< integer > fractureState( 1 );

  fractureState[0] = fields::contact::FractureState::Stick;
  array2d< real64 > traction( 1, 3 );
  array2d< real64 > jump( 1, 3 );
  array2d< real64 > djump( 1, 3 );

  bool isSimultaneous = m_isSimultaneous;
  real64 const slidingCheckTol = .05; //default
  real64 cos_phi = cos( m_phi * M_PI/180 );
  real64 cos_theta = cos( m_theta * M_PI/180 );
  real64 sin_phi = sin( m_phi * M_PI/180 );
  real64 sin_theta = sin( m_theta * M_PI/180 );

  real64 normalDispTolFac = m_normalDispTolFac;
  real64 normalTracTolFac = m_normalTracTolFac;
  real64 slidingTolFac = m_slidingTolFac;
  real64 iterPenNFac = m_iterPenNFac;
  real64 iterPenTFac = m_iterPenTFac;

  real64 area = m_area;
  real64 volumes[2]{m_volume[0], m_volume[1]};
  real64 shearModuli[2]{m_shearModulus[0], m_shearModulus[1]};
  real64 bulkModuli[2]{m_bulkModulus[0], m_bulkModulus[1]};

  // mimicAssembly only models a single homogeneous material for its two-hex
  // mock geometry -- average the two neighboring cells' moduli.
  array2d<real64> const D = mimicAssembly::buildD(
      0.5*(bulkModuli[0]+bulkModuli[1]),
      0.5*(shearModuli[0]+shearModuli[1]) );

  // Persistent contact state carried across rows. theta/phi are fixed
  // geometry inputs, so Kelastic only needs to be built once.
  mimicAssembly::ContactState contact{
      table( 0, NTRAC ),   // initial normal traction guess
      true,                                     // start unconstrained
      m_theta * M_PI/180.0,
      m_phi   * M_PI/180.0 };

  mimicAssembly::CRSMat const Kelastic = mimicAssembly::assemble( D, contact );

  real64 cumulativeImposedUx = 0.0;   // running sum of table(ei,NDJUMP)

  //TODO computeTolerance eleme to Elem
  typename FRICTION_TYPE::KernelWrapper const kernelWrapper = friction.createKernelUpdates();

  integer const numRows = m_table.size( 0 );
  // forAll< parallelDevicePolicy<> >( numRows,
  //                                   [&friction, &table, &kernelWrapper,
  //                                    &ghostRank,
  //                                    //  &normalDisplacementTol, &normalTractionTol, &slidingTol,
  //                                    &normalDispTolFac, &normalTracTolFac, &slidingTolFac,
  //                                    &iterPenNFac, &iterPenTFac,
  //                                    &area, &volumes, &bulkModuli, &shearModuli,
  //                                    &cos_phi, &sin_phi, &cos_theta, &sin_theta,
  //                                    &slidingCheckTol, &isSimultaneous,
  //                                    &jump, &djump,
  //                                    &fractureState, &traction ]
  //                                   GEOS_HOST_DEVICE ( integer const ei )
  
  for(integer ei=0; ei<numRows;ei+=20)
  {

    GEOS_LOG_RANK( "[debug] Table Evaluation" );

    jump[0][0] = table( ei, NJUMP );
    jump[0][1] = table( ei, SLIP0 );
    jump[0][2] = table( ei, SLIP1 );

    djump[0][0] = table( ei, NDJUMP );
    djump[0][1] = table( ei, DSLIP0 );
    djump[0][2] = table( ei, DSLIP1 );

    traction[0][0] = table( ei, NTRAC );
    traction[0][1] = table( ei, STRAC0 );
    traction[0][2] = table( ei, STRAC1 );

    kernelWrapper.updateFractureState( 0, 
                                       jump[0],
                                       traction[0],
                                       fractureState[0] );

    real64 normalDispTol{}, slidingTol{}, normalTractionTol{};
    real64 iterativePen[2]{};


    array2d< real64 > rotationMatrix( 3, 3 );
    rotationMatrix[0][0] = cos_phi*cos_theta;  rotationMatrix[0][1] = -sin_phi;  rotationMatrix[0][2] = cos_phi*sin_theta;
    rotationMatrix[1][0] = sin_phi*sin_theta;  rotationMatrix[1][1] = cos_phi;  rotationMatrix[1][2] = sin_theta*sin_phi;
    rotationMatrix[2][0] = -sin_theta;  rotationMatrix[2][1] = 0.;  rotationMatrix[2][2] = cos_theta;

    SolidMechanicsAugmentedLagrangianContact::computeTolerancePerFace( area,
                                                                       volumes,
                                                                       bulkModuli,
                                                                       shearModuli,
                                                                       rotationMatrix,
                                                                       normalDispTolFac,
                                                                       normalTracTolFac,
                                                                       slidingTolFac,
                                                                       iterPenNFac,
                                                                       iterPenTFac,
                                                                       normalDispTol,
                                                                       slidingTol,
                                                                       normalTractionTol,
                                                                       iterativePen );

    //small adapter
    array2d< real64 > const iterPen( 1, 5 );
    iterPen[0][0] = iterativePen[0];
    iterPen[0][1] = iterativePen[1];
    iterPen[0][2] = iterativePen[1]; iterPen[0][3] = iterativePen[1]; iterPen[0][4] = 0.;

    array1d< real64 > const a_normalDisplacementTol( 1 ); a_normalDisplacementTol[0]=normalDispTol;//normalDispTol should scale as 1/E
    array1d< real64 > const a_normalTractionTol( 1 ); a_normalTractionTol[0]=normalTractionTol;
    array1d< real64 > const a_slidingTol( 1 ); a_slidingTol[0]=slidingTol;


    

    auto [newTraction, condCov] = SolidMechanicsAugmentedLagrangianContact::updateTractionAndConstraintCheck( 1,
                                                                                                              friction,
                                                                                                              isSimultaneous,
                                                                                                              slidingCheckTol,
                                                                                                              a_normalDisplacementTol,
                                                                                                              a_normalTractionTol,
                                                                                                              a_slidingTol,
                                                                                                              iterPen,
                                                                                                              jump,
                                                                                                              djump,
                                                                                                              ghostRank,
                                                                                                              fractureState.toView(),
                                                                                                              traction.toView()
                                                                                                              );

    kernelWrapper.updateFractureState( 0, 
                                       jump[0],
                                       newTraction[0],
                                       fractureState[0] );




    cumulativeImposedUx += table( ei, DISP );

  // Bridges solveStep's Newton iterate (tNold, gN) to the actual configured
  // friction model, reusing the same call as the primary update above. Built
  // fresh each row: captures this row's tolerance/iterPen locals by reference.
  std::function< mimicAssembly::ContactState(double,double,double) > updateNormalTraction =
    [&]( double tNold, double gN, double /*epsNArg*/ ) -> mimicAssembly::ContactState
  {
    jump[0][0] = gN;         jump[0][1] = 0.0;                    jump[0][2] = 0.0;
    traction[0][0] = tNold;  traction[0][1] = table( ei, STRAC0 ); traction[0][2] = table( ei, STRAC1 );

    kernelWrapper.updateFractureState( 0, jump[0], traction[0], fractureState[0] );

    std::tie(newTraction, condCov) = SolidMechanicsAugmentedLagrangianContact::updateTractionAndConstraintCheck(
        1, friction, isSimultaneous, slidingCheckTol,
        a_normalDisplacementTol, a_normalTractionTol, a_slidingTol,
        iterPen, jump, djump, ghostRank, fractureState.toView(), traction.toView() );

    kernelWrapper.updateFractureState( 0, jump[0], newTraction[0], fractureState[0] );

    return mimicAssembly::ContactState{
        newTraction[0][0],
        fractureState[0] == fields::contact::FractureState::Open,
        contact.theta, contact.phi };
  };

  // epsN: reuse this row's normal iterative penalty (already folds in
  // area/volume/moduli scaling via computeTolerancePerFace) so the FEM
  // Newton loop's contact tangent stays consistent with the table-driven path.
  const auto stepResult =
      mimicAssembly::solveStep( cumulativeImposedUx, contact, Kelastic, iterativePen[0], updateNormalTraction );

  for(const auto& step : stepResult ){
  table( ei, FEMNTRAC )      = contact.tN;
  table( ei, FEMGN )         = step.gN;
  table( ei, FEMNEWTONITER ) = step.newtonIterations;
  table( ei, FEMCONVERGED )  = step.converged ? 1.0 : 0.0;
      
  //old cols
  table( ei, CC ) = condCov[0];
  table( ei, FS ) = fractureState[0];
  table( ei, NEWTRAC ) = newTraction[0][0];
  table( ei, SNEWTRAC0 ) = newTraction[0][1];
  table( ei, SNEWTRAC1 ) = newTraction[0][2];


  table( ei, NTOL ) = normalDispTol;
  table( ei, TTOL ) = slidingTol;
  table( ei, NTRACTOL ) = normalTractionTol;

  table( ei, ITERPEN0 ) = iterPen[0][0];
  table( ei, ITERPEN1 ) = iterPen[0][1];

  }
  // } ); //TODO restore once forAll<>
}
}//function
}//namespace


#endif //GEOS_FRICTIONDRIVERRUNTEST_HPP_
