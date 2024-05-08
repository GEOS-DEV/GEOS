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
 * @file CompositionalMultiphaseFVMApplyFaceDirichletBC.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseFVMKernels.hpp"

namespace geos
{

using namespace constitutive;

void CompositionalMultiphaseFVM::applyFaceDirichletBC( real64 const time_n,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateFaceDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getDataContext() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    applyFieldValue< FaceManager >( time_n, dt, mesh, getFaceBoundaryConditionMessage(),
                                    fields::flow::pressure::key(), fields::flow::facePressure::key() );
    // Take BCs defined for "globalCompFraction" field and apply values to "faceGlobalCompFraction"
    applyFieldValue< FaceManager >( time_n, dt, mesh, getFaceBoundaryConditionMessage(),
                                    fields::flow::globalCompFraction::key(), fields::flow::faceGlobalCompFraction::key() );
    // Take BCs defined for "temperature" field and apply values to "faceTemperature"
    applyFieldValue< FaceManager >( time_n, dt, mesh, getFaceBoundaryConditionMessage(),
                                    fields::flow::temperature::key(), fields::flow::faceTemperature::key() );

    // Then launch the face Dirichlet kernel
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    fields::flow::pressure::key(), // we have required that pressure is always present
                                    [&] ( FieldSpecificationBase const &,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const &,
                                          FaceManager &,
                                          string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( stencil.size() == 0 )
      {
        return;
      }

      // TODO: same issue as in the single-phase case
      //       currently we just use model from the first cell in this stencil
      //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
      //       Can we just use cell properties for an approximate flux computation?
      //       Then we can forget about capturing the fluid model.
      localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
      localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
      ElementSubRegionBase & subRegion = elemManager.getRegion( er ).getSubRegion( esr );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & multiFluidBase = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

      string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        //todo (jafranc) extend upwindScheme name if satisfied in isothermalCase
        thermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        isothermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }

    } );
  } );
}

}// namespace geos
