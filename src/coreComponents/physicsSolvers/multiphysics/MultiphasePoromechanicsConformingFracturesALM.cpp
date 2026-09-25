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
 * @file MultiphasePoromechanicsConformingFracturesALM.cpp
 */

#include "MultiphasePoromechanicsConformingFracturesALM.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::MultiphasePoromechanicsConformingFracturesALM( const string & name,
                                                                                                             Group * const parent )
  : Base( name, parent )
{}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setSparsityPattern( DomainPartition & domain,
                                                                                       DofManager & dofManager,
                                                                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                       SparsityPattern< globalIndex > & pattern )
{
  GEOS_MARK_FUNCTION;

  // Recompute fracture face/element geometry and rebuild the ALM contact solver's internal lists.
  // These must happen before assembling the contact-dependent pattern; setSparsityPattern() is the
  // only place this coupled solver hooks into for that, since it does not route through the contact
  // sub-solver's own setupSystem().
  this->solidMechanicsSolver()->updateFractureGeometry( domain );

  Base::setSparsityPattern( domain, dofManager, localMatrix, pattern );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( string const & GEOS_UNUSED_PARAM( meshName ),
                                            MeshLevel const & GEOS_UNUSED_PARAM( mesh ),
                                            string_array const & GEOS_UNUSED_PARAM( regionNames ),
                                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                            CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                            arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  // Only reached through the monolithic assembleSystem(), which is exercised by FullyImplicit
  // coupling. The compositional-flow ALM force/pressure coupling kernels are not implemented yet;
  // this solver is sequential-only for now (see postInputInitialization/couplingType checks upstream).
  GEOS_ERROR( GEOS_FMT( "{}: FullyImplicit coupling is not supported by {}", this->getName(), this->getCatalogName() ) );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( string const & GEOS_UNUSED_PARAM( meshName ),
                                                    MeshLevel const & GEOS_UNUSED_PARAM( mesh ),
                                                    string_array const & GEOS_UNUSED_PARAM( regionNames ),
                                                    DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                    CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  // See assembleForceResidualDerivativeWrtPressure: FullyImplicit-only path, not implemented for
  // compositional flow yet.
  GEOS_ERROR( GEOS_FMT( "{}: FullyImplicit coupling is not supported by {}", this->getName(), this->getCatalogName() ) );
}

template class MultiphasePoromechanicsConformingFracturesALM<>;
template class MultiphasePoromechanicsConformingFracturesALM< CompositionalMultiphaseReservoirAndWells<> >;

namespace
{
typedef MultiphasePoromechanicsConformingFracturesALM< CompositionalMultiphaseReservoirAndWells<> > MultiphaseReservoirPoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphaseReservoirPoromechanicsConformingFracturesALM, string const &, Group * const )
typedef MultiphasePoromechanicsConformingFracturesALM<> MultiphasePoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphasePoromechanicsConformingFracturesALM, string const &, Group * const )
}

} /* namespace geos */
