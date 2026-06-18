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
 * @file SinglePhasePoromechanicsConformingFracturesALM.cpp
 */

#include "SinglePhasePoromechanicsConformingFracturesALM.hpp"

#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::SinglePhasePoromechanicsConformingFracturesALM( const string & name,
                                                                                                               Group * const parent )
  : Base( name, parent )
{}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::postInputInitialization()
{
  Base::postInputInitialization();

  GEOS_WARNING_IF( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit,
                   "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );

  this->getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                                   DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  Base::setupCoupling( domain, dofManager );

  dofManager.addCoupling( this->getFlowDofKey(),
                          fields::contact::traction::key(),
                          DofManager::Connector::Elem );

}


template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupSystem( DomainPartition & domain,
                                                                                 DofManager & dofManager,
                                                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                 ParallelVector & rhs,
                                                                                 ParallelVector & solution,
                                                                                 bool const setSparsity )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                    real64 const dt,
                                                                                    DomainPartition & domain,
                                                                                    DofManager const & dofManager,
                                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                    arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleElementBasedContributions( real64 const time_n,
                                                                                                       real64 const dt,
                                                                                                       DomainPartition & domain,
                                                                                                       DofManager const & dofManager,
                                                                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                       arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleCouplingTerms( real64 const time_n,
                                                                                           real64 const dt,
                                                                                           DomainPartition const & domain,
                                                                                           DofManager const & dofManager,
                                                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, localRhs, time_n, dt );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  Base::updateState( domain );
  this->solidMechanicsSolver()->updateState( domain );

  this->flowSolver()->prepareStencilWeights( domain );
  updateHydraulicApertureAndFracturePermeability( domain );
  this->flowSolver()->updateStencilWeights( domain );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
setUpDflux_dApertureMatrix( DomainPartition & domain,
                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                            CRSMatrix< real64, globalIndex > & localMatrix )
{
  GEOS_UNUSED_VAR( domain, localMatrix );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                DofManager const & dofManager,
                                arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, rowLengths );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( domain, dofManager, pattern );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                            MeshLevel const & mesh,
                                            arrayView1d< string const > const & regionNames,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    arrayView1d< string const > const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( mesh, regionNames, dofManager, localMatrix );

  if( this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential )
  {
    GEOS_ERROR( "SinglePhasePoromechanicsConformingFracturesALM does not support FullyImplicit coupling type." );
  }

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::mapSolutionBetweenSolvers( DomainPartition & domain,
                                                                                                integer const solverType )
{
  GEOS_MARK_FUNCTION;

  if( solverType == static_cast< integer >( Base::SolverType::SolidMechanics )
      && !this->m_performStressInitialization )
  {
    this->flowSolver()->prepareStencilWeights( domain );
    updateHydraulicApertureAndFracturePermeability( domain );
    this->flowSolver()->updateStencilWeights( domain );
  }

  Base::mapSolutionBetweenSolvers( domain, solverType );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateHydraulicApertureAndFracturePermeability( DomainPartition & domain )
{
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const dispJump = subRegion.getField< fields::contact::dispJump >();
      arrayView1d< real64 const > const area = subRegion.getElementArea();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView2d< real64 const > const fractureTraction = subRegion.getField< fields::contact::traction >();
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

      arrayView1d< real64 > const aperture = subRegion.getElementAperture();
      arrayView1d< real64 > const hydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
      arrayView1d< real64 > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();
      arrayView1d< integer > const & fractureState = subRegion.getField< fields::contact::fractureState >();

      string const porousSolidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      string const & hydraulicApertureRelationName = subRegion.getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString() );
      HydraulicApertureBase const & hydraulicApertureModel =
        this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

      constitutiveUpdatePassThru( hydraulicApertureModel, [&] ( auto & castedHydraulicAperture )
      {
        using HydraulicApertureType = TYPEOFREF( castedHydraulicAperture );
        typename HydraulicApertureType::KernelWrapper hydraulicApertureWrapper = castedHydraulicAperture.createKernelWrapper();

        ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
        {
          typename TYPEOFREF( castedPorousSolid )::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();

          poromechanicsFracturesKernels::StateUpdateKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                                               porousMaterialWrapper,
                                                                                               hydraulicApertureWrapper,
                                                                                               dispJump,
                                                                                               pressure,
                                                                                               area,
                                                                                               volume,
                                                                                               deltaVolume,
                                                                                               aperture,
                                                                                               oldHydraulicAperture,
                                                                                               hydraulicAperture,
                                                                                               fractureTraction,
                                                                                               fractureState );
        } );
      } );
    } );
  } );

}


template class SinglePhasePoromechanicsConformingFracturesALM<>;
template class SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> >;

namespace
{
typedef SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanicsConformingFracturesALM, string const &, Group * const )
typedef SinglePhasePoromechanicsConformingFracturesALM<> SinglePhasePoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsConformingFracturesALM, string const &, Group * const )
}

} /* namespace geos */
