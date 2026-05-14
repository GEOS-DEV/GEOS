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
 * @file CompositionalMultiphaseFVM_applyFaceDirichletBC.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalDirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/DirichletFluxComputeZFormulationKernel.hpp"

namespace geos
{

using namespace constitutive;
using namespace fields;

namespace
{
constexpr const char* faceBcLogMessage =
  "CompositionalMultiphaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target faces (including ghost faces) is {}."
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

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
    GEOS_ERROR_IF( !bcConsistent, "inconsistent boundary conditions", getDataContext() );
  }

  using namespace isothermalCompositionalMultiphaseFVMKernels;

  // for now, we neglect capillary pressure in the kernel
  BitFlags< KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( KernelFlags::TotalMassEquation );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::pressure::key(), flow::facePressure::key() );
    // Take BCs defined for "globalCompFraction" field and apply values to "faceGlobalCompFraction"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::globalCompFraction::key(), flow::faceGlobalCompFraction::key() );
    // Take BCs defined for "temperature" field and apply values to "faceTemperature"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::temperature::key(), flow::faceTemperature::key() );

    // Then launch the face Dirichlet kernel
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    flow::pressure::key(), // we have required that pressure is always present
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

      if( m_formulationType == CompositionalMultiphaseFormulationType::OverallComposition )
      {
        isothermalCompositionalMultiphaseFVMKernels::
          DirichletFluxComputeZFormulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     kernelFlags,
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
        if( m_isThermal )
        {
          //todo (jafranc) extend upwindScheme name if satisfied in isothermalCase
          thermalCompositionalMultiphaseFVMKernels::
            DirichletFluxComputeKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       kernelFlags,
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
            DirichletFluxComputeKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       kernelFlags,
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
      }

    } );
  } );
}

}// namespace geos
