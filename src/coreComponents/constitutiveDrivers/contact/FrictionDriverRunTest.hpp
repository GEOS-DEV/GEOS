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

#include <conduit.hpp>

namespace geos
{

template< typename FRICTION_TYPE >
void
FrictionDriver::runTest( FRICTION_TYPE & friction,
                         const arrayView2d< real64 > & table )
{

  array1d< integer > const ghostRank( 1 ); ghostRank[0] = -1;
  array1d< real64 > const normalDisplacementTol( 1 ); normalDisplacementTol[0]=10;
  array1d< real64 > const normalTractionTol( 1 ); normalTractionTol[0]=1.;
  array1d< real64 > const slidingTol( 1 ); slidingTol[0]=.1;//normalDispTol should scale as 1/E

  real64 kt = 1.;
  array2d< real64 > const iterPen( 1, 5 );
  iterPen[0][0] = 100*kt;
  iterPen[0][1] = kt;
  iterPen[0][2] = kt; iterPen[0][3] = kt; iterPen[0][4] = 0.;
  array1d< integer > fractureState( 1 );

  fractureState[0] = fields::contact::FractureState::Stick;
  array2d< real64 > traction( 1, 3 );
  array2d< real64 > jump( 1, 3 );
  array2d< real64 > djump( 1, 3 );

  //TODO computeTolerance eleme to Elem

  integer const numRows = m_table.size( 0 );
  forAll< parallelDevicePolicy<> >( numRows,
                                    [&friction, &table,
                                     &ghostRank, &normalDisplacementTol, &normalTractionTol, &slidingTol, &iterPen,
                                     &jump, &djump,
                                     &fractureState, &traction ]
                                    GEOS_HOST_DEVICE ( integer const ei )
  {

    GEOS_LOG_RANK( "[debug] Table Evaluation" );
    GEOS_LOG_RANK( GEOS_FMT( " jump [{}x{}]", jump.size( 0 ), jump.size( 1 )));

    jump[0][0] = table( ei, NJUMP );
    jump[0][1] = table( ei, SLIP0 );
    jump[0][2] = table( ei, SLIP1 );

    djump[0][0] = table( ei, NDJUMP );
    djump[0][1] = table( ei, DSLIP0 );
    djump[0][2] = table( ei, DSLIP1 );

    traction[0][0] = table( ei, NTRAC );
    traction[0][1] = table( ei, STRAC0 );
    traction[0][2] = table( ei, STRAC1 );


    SolidMechanicsAugmentedLagrangianContact::updateTractionAndConstraintCheck( 1,
                                                                                friction,
                                                                                true,//simultaneous
                                                                                0.05,//slidingCheckTolerance
                                                                                normalDisplacementTol,
                                                                                normalTractionTol,
                                                                                slidingTol,
                                                                                iterPen,
                                                                                jump,
                                                                                djump,
                                                                                ghostRank,
                                                                                fractureState.toView(),
                                                                                traction.toView()
                                                                                );
    table( ei, FS ) = fractureState[0];
    table( ei, NEWTRAC ) = traction[0][0];
    table( ei, SNEWTRAC0 ) = traction[0][1];
    table( ei, SNEWTRAC1 ) = traction[0][2];

  } );
}

}

#endif //GEOS_FRICTIONDRIVERRUNTEST_HPP_
