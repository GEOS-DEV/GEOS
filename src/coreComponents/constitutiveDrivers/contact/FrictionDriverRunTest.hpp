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
#include "constitutive/solid/SolidFields.hpp"

#include "mesh/DomainPartition.hpp"
#include <conduit.hpp>

namespace geos
{

template< typename FRICTION_TYPE, typename CONTACT_SOLVER >
void
FrictionDriver::runTest( FRICTION_TYPE & friction,
                         CONTACT_SOLVER & contact,
                         const arrayView2d< real64 > & table )
{
  // Create kernel wrapper and trigger solver configuration update.
  conduit::Node dummyRoot;
  DomainPartition dummyDomain( "FrictionDriverRunTestDomain", dummyRoot );
  // CONTACT_SOLVER const solver( "FrictionDriverRunTestSolver", &dummyDomain );
  // solver.updateConfiguration( dummyDomain, 0 );
  contact.updateConfiguration( dummyDomain, 0 );

  typename FRICTION_TYPE::KernelWrapper const kernelWrapper = friction.createKernelUpdates();

  integer const numRows = m_table.size( 0 );
  forAll< parallelDevicePolicy<> >( numRows,
                                    [ kernelWrapper, table ]
                                    GEOS_HOST_DEVICE ( integer const ei )
  {
    stackArray1d< real64, 3 > jump( 3 );
    stackArray1d< real64, 3 > traction( 3 );

    jump[0] = table( ei, NJUMP );
    jump[1] = table( ei, SLIP0 );
    jump[2] = table( ei, SLIP1 );

    traction[0] = table( ei, NTRAC );
    traction[1] = table( ei, STRAC0 );
    traction[2] = table( ei, STRAC1 );

    integer fracture_state = fields::contact::FractureState::Stick;
    kernelWrapper.updateFractureState( jump.toSliceConst(),
                                       traction.toSliceConst(),
                                       fracture_state );

    table( ei, FS ) = fracture_state;
  } );
}

}

#endif //GEOS_FRICTIONDRIVERRUNTEST_HPP_
