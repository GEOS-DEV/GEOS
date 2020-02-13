/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseFVM.cpp
 */

#include "SinglePhaseFVM.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "common/TimingMacros.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseKernels;

SinglePhaseFVM::SinglePhaseFVM( const std::string& name,
                                Group * const parent ):
  SinglePhaseBase(name, parent)
{
  m_numDofPerCell = 1;
}


void SinglePhaseFVM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::pressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connectivity::Face );
}
  

real64 SinglePhaseFVM::CalculateResidualNorm( DomainPartition const * const domain,
                                              DofManager const & dofManager,
                                              ParallelVector const & rhs )
{
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm[3] = { 0.0, 0.0, 0.0 };
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64 const> const & refPoro        = m_porosityRef[er][esr];
    arrayView1d<real64 const> const & volume         = m_volume[er][esr];
    arrayView1d<real64 const> const & densOld        = m_densityOld[er][esr];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {
      if (elemGhostRank[a] < 0)
      {
        localIndex const lid = rhs.getLocalRowID( dofNumber[a] );
        real64 const val = localResidual[lid];
        localResidualNorm[0] += val * val;
        localResidualNorm[1] += refPoro[a] * densOld[a] * volume[a];
        localResidualNorm[2] += 1;
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm[3] = {0,0,0};
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         3,
                         MPI_SUM,
                         MPI_COMM_GEOSX);


  real64 const residual = sqrt(globalResidualNorm[0]) / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid ) = (%4.2e) ; ",
             residual);
    std::cout<<output;
  }


  return residual;
}

  
void SinglePhaseFVM::ApplySystemSolution( DofManager const & dofManager,
                                          ParallelVector const & solution,
                                          real64 const scalingFactor,
                                          DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  dofManager.addVectorToField( solution,
                               viewKeyStruct::pressureString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor );

  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );
}

void SinglePhaseFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                        real64 const dt,
                                        DomainPartition const * const domain,
                                        DofManager const * const dofManager,
                                        ParallelMatrix * const matrix,
                                        ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager=  mesh->getElemManager();

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > dofNumberAccessor =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofKey );

  FluxKernel::ElementView< arrayView1d<globalIndex const> > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & gravCoef   = m_gravCoef.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens        = m_density.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres = m_dDens_dPres.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & mob         = m_mobility.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres  = m_dMobility_dPres.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture0  = m_elementAperture0.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture  = m_effectiveAperture.toViewConst();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  FluxKernel::ElementView < arrayView1d<real64 const> > const & separationCoeff  = m_elementSeparationCoefficient.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & dseparationCoeff_dAper  = m_element_dSeparationCoefficient_dAperture.toViewConst();
#endif
  localIndex const fluidIndex = m_fluidIndex;


  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

//    typedef TYPEOFREF( stencil ) STENCIL_TYPE;

    FluxKernel::Launch( stencil,
                        dt,
                        fluidIndex,
                        dofNumber,
                        pres,
                        dPres,
                        gravCoef,
                        dens,
                        dDens_dPres,
                        mob,
                        dMob_dPres,
                        aperture0,
                        aperture,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                        separationCoeff,
                        dseparationCoeff_dAper,
#endif
                        matrix,
                        rhs,
                        *m_derivativeFluxResidual_dAperture );
  });
}

void
SinglePhaseFVM::ApplyBoundaryConditions( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition * const domain,
                                         DofManager const & dofManager,
                                         ParallelMatrix & matrix,
                                         ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  // call the BoundaryConditionManager::ApplyField function that will check to see
  // if the boundary condition should be applied to this subregion
  fsManager.Apply( time_n + dt, domain, "ElementRegions", FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        set<localIndex> const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d< integer const > const &
    ghostRank = subRegion->getReference<array1d<integer> >( ObjectManagerBase::viewKeyStruct::ghostRankString);

    set< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert(a);
      }
    }

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( localSet,
                                                                            time_n + dt,
                                                                            dt,
                                                                            subRegion,
                                                                            dofNumber,
                                                                            1,
                                                                            matrix,
                                                                            rhs,
                                                                            [&]( localIndex const GEOSX_UNUSED_ARG( a ) ) -> real64
    {
      return 0;
    } );

  } );


  fsManager.Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        set<localIndex> const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    //for now assume all the non-flux boundary conditions are Dirichlet type BC.

    arrayView1d<real64 const> const &
    pres = subRegion->getReference<array1d<real64> >( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const &
    dPres = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                              time_n + dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              1,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const a ) -> real64
    {
      return pres[a] + dPres[a];
    });
  });


  ApplyFaceDirichletBC_implicit( time_n, dt, &dofManager, domain, &matrix, &rhs );

  matrix.close();
  rhs.close();

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After SinglePhaseFVM::ApplyBoundaryConditions" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOSX_LOG_RANK_0( "After SinglePhaseFVM::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void SinglePhaseFVM::ApplyFaceDirichletBC_implicit( real64 const time_n,
                                                    real64 const dt,
                                                    DofManager const * const dofManager,
                                                    DomainPartition * const domain,
                                                    ParallelMatrix * const matrix,
                                                    ParallelVector * const rhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  arrayView2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();

  ConstitutiveManager * const constitutiveManager =
    domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  NumericalMethodsManager * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager * const fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > dofNumberAccessor =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofKey );

  FluxKernel::ElementView< arrayView1d<globalIndex const> > const & dofNumber = dofNumberAccessor.toViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & dPres       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & gravCoef    = m_gravCoef;
  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & dens        = m_density;
  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & mob         = m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & dMob_dPres  = m_dMobility_dPres;

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // use ArrayView to make capture by value easy in lambdas
  arrayView1d<real64 const> const & presFace      = faceManager->getReference< array1d<real64> >( viewKeyStruct::boundaryFacePressureString );
  arrayView2d<real64>       const & densFace      = faceManager->getReference< array2d<real64> >( viewKeyStruct::boundaryFaceDensityString );
  arrayView2d<real64>       const & viscFace      = faceManager->getReference< array2d<real64> >( viewKeyStruct::boundaryFaceViscosityString );
  arrayView1d<real64>       const & mobFace       = faceManager->getReference< array1d<real64> >( viewKeyStruct::boundaryFaceMobilityString );
  arrayView1d<real64 const> const & gravCoefFace  = faceManager->getReference< array1d<real64> >( viewKeyStruct::gravityCoefString );

  dataRepository::Group const * sets = faceManager->sets();

  // first, evaluate BC to get primary field values (pressure)
//  fsManager->ApplyField(faceManager, viewKeyStruct::boundaryFacePressure, time + dt);
  fsManager.Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::boundaryFacePressureString,
                    [&] ( FieldSpecificationBase const * const fs,
                          string const &,
                          set<localIndex> const & targetSet,
                          Group * const targetGroup,
                          string const fieldName )
  {
    fs->ApplyFieldValue<FieldSpecificationEqual>(targetSet,time_n + dt, targetGroup, fieldName);
  });


  // call constitutive models to get dependent quantities needed for flux (density, viscosity)
  fsManager.Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::boundaryFacePressureString,
                    [&] ( FieldSpecificationBase const * GEOSX_UNUSED_ARG( bc ),
                          string const &,
                          set<localIndex> const & targetSet,
                          Group * const,
                          string const & )
  {
    for (auto kf : targetSet)
    {
      // since we don't have models on faces yet, we take them from an adjacent cell
      integer ke;
      for (ke = 0; ke < 2; ++ke)
      {
        if (elemRegionList[kf][ke] >= 0 && regionFilter.contains(elemRegionList[kf][ke]))
        {
          break;
        }
      }
      GEOSX_ERROR_IF( ke > 1, "Face not adjacent to target regions: " << kf );
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];

      real64 dummy; // don't need derivatives on faces

      SingleFluidBase * fluid = constitutiveRelations[er][esr][m_fluidIndex]->group_cast<SingleFluidBase *>();
      fluid->Compute( presFace[kf], densFace[kf][0], dummy, viscFace[kf][0], dummy );
    }

    MobilityKernel::Launch( targetSet, densFace, viscFace, mobFace );
  });

  // *** assembly loop ***

  constexpr localIndex numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  real64 densWeight[numElems] = { 0.5, 0.5 };

  fsManager.Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::boundaryFacePressureString,
                    [&] ( FieldSpecificationBase const * GEOSX_UNUSED_ARG( bc ),
                          string const & setName,
                          set<localIndex> const &,
                          Group * const,
                          string const & )
  {
    if (!sets->hasWrapper(setName) || !fluxApprox->hasBoundaryStencil(setName))
      return;

    FluxApproximationBase::BoundaryStencil const & stencil = fluxApprox->getBoundaryStencil(setName);
    ArrayOfArraysView<FluxApproximationBase::BoundaryStencil::Entry const, true> const & connections = stencil.getConnections();

    forall_in_range<serialPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray(iconn);

      stackArray1d<globalIndex, maxStencilSize> dofColIndices( stencilSize );

      stackArray1d<real64, numElems> mobility( numElems );
      stackArray1d<real64, numElems> dMobility_dP( numElems );
      stackArray1d<real64, maxStencilSize> dDensMean_dP( stencilSize );
      stackArray1d<real64, maxStencilSize> dFlux_dP( stencilSize );
      stackArray1d<real64, maxStencilSize> localFluxJacobian( stencilSize );

      // clear working arrays
      dDensMean_dP = 0.0;

      // calculate quantities on primary connected points
      real64 densMean = 0.0;
      globalIndex eqnRowIndex = -1;
      localIndex cell_order = -1;

      for (localIndex i = 0; i < numElems; ++i)
      {
        PointDescriptor const & point = connections(iconn, i).index;

        real64 density = 0, dDens_dP = 0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er  = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei  = point.cellIndex.index;

            eqnRowIndex = dofNumber[er][esr][ei];

            density  = dens[er][esr][m_fluidIndex][ei][0];
            dDens_dP = dDens_dPres[er][esr][m_fluidIndex][ei][0];

            mobility[i]     = mob[er][esr][ei];
            dMobility_dP[i] = dMob_dPres[er][esr][ei];

            cell_order = i; // mark position of the cell in connection for sign consistency later
            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            density  = densFace[point.faceIndex][0];
            dDens_dP = 0.0;

            mobility[i]     = mobFace[point.faceIndex];
            dMobility_dP[i] = 0.0;
            break;
          }
          default:
            GEOSX_ERROR("Unsupported point type in stencil");
        }

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      dofColIndices = -1;
      for (localIndex i = 0; i < stencilSize; ++i)
      {
        FluxApproximationBase::BoundaryStencil::Entry const & entry = connections(iconn, i);
        PointDescriptor const & point = entry.index;

        real64 pressure = 0.0, gravD = 0.0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei = point.cellIndex.index;

            dofColIndices[i] = dofNumber[er][esr][ei];
            pressure = pres[er][esr][ei] + dPres[er][esr][ei];
            gravD = gravCoef[er][esr][ei];

            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            localIndex const kf = point.faceIndex;

            pressure = presFace[kf];
            gravD = gravCoefFace[kf];

            break;
          }
          default:
          GEOSX_ERROR("Unsupported point type in stencil");
        }

        real64 const gravTerm = densMean * gravD;
        real64 const dGrav_dP = dDensMean_dP[i] * gravD;

        potDif += entry.weight * (pressure + gravTerm);
        dFlux_dP[i] = entry.weight * (1.0 + dGrav_dP);
      }

      // upwinding of fluid properties (make this an option?)
      localIndex const k_up = (potDif >= 0) ? 0 : 1;

      // compute the final flux and derivatives
      real64 const flux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
        dFlux_dP[ke] *= mobility[k_up];
      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      integer sign = (cell_order == 0 ? 1 : -1);
      real64 const localFlux =  dt * flux * sign;

      integer counter = 0;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        // compress arrays, skipping face derivatives
        if (dofColIndices[ke] >= 0)
        {
          dofColIndices[counter] = dofColIndices[ke];
          localFluxJacobian[counter] = dt * dFlux_dP[ke] * sign;
          ++counter;
        }
      }

      // Add to global residual/jacobian
      matrix->add( eqnRowIndex, dofColIndices.data(), localFluxJacobian.data(), counter );
      rhs->add( eqnRowIndex, localFlux );
    } );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFVM, std::string const &, Group * const )
} /* namespace geosx */
