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
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                                  DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager );

  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );

}


template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupSystem( DomainPartition & domain,
                                                                                DofManager & dofManager,
                                                                                CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                ParallelVector & rhs,
                                                                                ParallelVector & solution,
                                                                                bool const setSparsity )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );

}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                   real64 const dt,
                                                                                   DomainPartition & domain,
                                                                                   DofManager const & dofManager,
                                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                   arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs );

  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleElementBasedContributions( real64 const time_n,
                                                                                                      real64 const dt,
                                                                                                      DomainPartition & domain,
                                                                                                      DofManager const & dofManager,
                                                                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                      arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs );

  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );

}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleCouplingTerms( real64 const time_n,
                                                                                          real64 const dt,
                                                                                          DomainPartition const & domain,
                                                                                          DofManager const & dofManager,
                                                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, localRhs, time_n, dt );

  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );
  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
}


template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                DofManager const & dofManager,
                                arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, rowLengths );
  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );

}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, pattern );
  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                            MeshLevel const & mesh,
                                            arrayView1d< string const > const & regionNames,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    arrayView1d< string const > const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( mesh, regionNames, dofManager, localMatrix );
  GEOS_ERROR( "MultiphasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );

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
