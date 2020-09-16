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
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;
using namespace SinglePhaseFVMKernels;

template< typename BASE >
SinglePhaseFVM< BASE >::SinglePhaseFVM( const std::string & name,
                                        Group * const parent ):
  BASE( name, parent )
{
  m_numDofPerCell = 1;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::SetupDofs( DomainPartition const * const domain,
                                        DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       targetRegionNames() );

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::pressureString, fluxApprox );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::SetupSystem( DomainPartition * const domain,
                                          DofManager & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  BASE::SetupSystem( domain,
                     dofManager,
                     matrix,
                     rhs,
                     solution );

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  std::unique_ptr< CRSMatrix< real64, localIndex > > &
  derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();
  {

    localIndex numRows = 0;
    forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & elementSubRegion )
    {
      numRows += elementSubRegion.size();
    } );

    derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );

    derivativeFluxResidual_dAperture->reserveNonZeros( matrix.numLocalNonzeros() );
    localIndex maxRowSize = -1;
    for( localIndex row=0; row<matrix.numLocalRows(); ++row )
    {
      localIndex const rowSize = matrix.localRowLength( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
    }
    for( localIndex row=matrix.numLocalRows(); row<numRows; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row,
                                                         maxRowSize );
    }

  }


  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  NumericalMethodsManager const *
    numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const *
    fvManager = numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( this->getDiscretization() );


  fluxApprox->forStencils< FaceElementStencil >( [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        for( localIndex k1=0; k1<numFluxElems; ++k1 )
        {
          derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
        }
      }
    }
  } );
}


template< typename BASE >
real64 SinglePhaseFVM< BASE >::CalculateResidualNorm( DomainPartition const * const domain,
                                                      DofManager const & dofManager,
                                                      ParallelVector const & rhs )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm[3] = { 0.0, 0.0, 0.0 };
  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d< real64 const > const & refPoro        = m_porosityRef[er][esr];
    arrayView1d< real64 const > const & volume         = m_volume[er][esr];
    arrayView1d< real64 const > const & densOld        = m_densityOld[er][esr];

    localIndex const subRegionSize = subRegion.size();
    for( localIndex a = 0; a < subRegionSize; ++a )
    {
      if( elemGhostRank[a] < 0 )
      {
        localIndex const lid = rhs.getLocalRowID( dofNumber[a] );
        real64 const val = localResidual[lid];
        localResidualNorm[0] += val * val;
        localResidualNorm[1] += refPoro[a] * densOld[a] * volume[a];
        localResidualNorm[2] += 1;
      }
    }
  } );

  // compute global residual norm
  real64 globalResidualNorm[3] = {0, 0, 0};
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         3,
                         MPI_SUM,
                         MPI_COMM_GEOSX );


  real64 const residual = sqrt( globalResidualNorm[0] ) / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid ) = (%4.2e) ; ",
             residual );
    std::cout<<output;
  }


  return residual;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::ApplySystemSolution( DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::pressureString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain->getNeighbors() );

  forTargetSubRegions( mesh, [&] ( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    this->UpdateState( subRegion, targetIndex );
  } );

}

template< typename BASE >
void SinglePhaseFVM< BASE >::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const dt,
                                                DomainPartition const * const domain,
                                                DofManager const * const dofManager,
                                                ParallelMatrix * const matrix,
                                                ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_derivativeFluxResidual_dAperture==nullptr )
  {
    m_derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( matrix->numLocalRows(),
                                                                                              matrix->numLocalCols() );
  }
  m_derivativeFluxResidual_dAperture->setValues( 0.0 );

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager=  mesh->getElemManager();

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d< globalIndex >,
                                                          arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementView< arrayView1d< globalIndex const > > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementView< arrayView1d< real64 const > > const & dPres       = m_deltaPressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & pres        = m_pressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & gravCoef    = m_gravCoef.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dens        = m_density.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dDens_dPres = m_dDens_dPres.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & mob         = m_mobility.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & dMob_dPres  = m_dMobility_dPres.toViewConst();

  FluxKernel::ElementView< arrayView1d< real64 const > > const & aperture0  = m_elementAperture0.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & aperture  = m_effectiveAperture.toViewConst();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  FluxKernel::ElementView< arrayView1d< real64 const > > const & separationCoeff  = m_elementSeparationCoefficient.toViewConst();

  FluxKernel::ElementView< arrayView1d< real64 const > > const & dseparationCoeff_dAper  = m_element_dSeparationCoefficient_dAperture.toViewConst();
#endif

  FluxKernel::ElementView< arrayView1d< R1Tensor const > > const & transTMultiplier  = m_transTMultiplier.toViewConst();

  fluxApprox->forAllStencils( [&]( auto const & stencil )
  {

//    typedef TYPEOFREF( stencil ) STENCIL_TYPE;

    FluxKernel::Launch( stencil,
                        dt,
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
                        transTMultiplier,
                        this->gravityVector(),
                        this->m_meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                        separationCoeff,
                        dseparationCoeff_dAper,
#endif
                        matrix,
                        rhs,
                        m_derivativeFluxResidual_dAperture->toView(),
			domain);
  } );
}

template< typename BASE >
void
SinglePhaseFVM< BASE >::ApplyBoundaryConditions( real64 const time_n,
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
  //TJ: add the source term (flux at elmt center)
  fsManager.Apply( time_n + dt, domain, "ElementRegions", FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const &
    ghostRank = subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    SortedArray< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert( a );
      }
    }

    fs->ApplyBoundaryConditionToSystem< FieldSpecificationAdd, LAInterface >( localSet.toViewConst(),
                                                                              time_n + dt,
                                                                              dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              1,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const GEOSX_UNUSED_PARAM( a ) ) -> real64
    {
      return 0;
    } );

  } );


  fsManager.Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    //for now assume all the non-flux boundary conditions are Dirichlet type BC.

    arrayView1d< real64 const > const &
    pres = subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );

    arrayView1d< real64 const > const &
    dPres = subRegion->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem< FieldSpecificationEqual, LAInterface >( lset,
                                                                                time_n + dt,
                                                                                subRegion,
                                                                                dofNumber,
                                                                                1,
                                                                                matrix,
                                                                                rhs,
                                                                                [&]( localIndex const a ) -> real64
    {
      return pres[a] + dPres[a];
    } );
  } );


  ApplyFaceDirichletBC_implicit( time_n, dt, &dofManager, domain, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After SinglePhaseFVM<BASE>::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string const filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string const filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After SinglePhaseFVM<BASE>::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

template< typename BASE >
void SinglePhaseFVM< BASE >::ApplyFaceDirichletBC_implicit( real64 const time_n,
                                                            real64 const dt,
                                                            DofManager const * const dofManager,
                                                            DomainPartition * const domain,
                                                            ParallelMatrix * const matrix,
                                                            ParallelVector * const rhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *mesh.getElemManager();
  FaceManager & faceManager = *mesh.getFaceManager();

  arrayView2d< localIndex > const & elemRegionList     = faceManager.elementRegionList();
  arrayView2d< localIndex > const & elemSubRegionList  = faceManager.elementSubRegionList();

  ConstitutiveManager * const constitutiveManager =
    domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  NumericalMethodsManager * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager * const fvManager = numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  // make a list of region indices to be included
  map< localIndex, localIndex > regionFluidMap;
  forTargetRegionsComplete( mesh, [&]( localIndex const targetIndex, localIndex const er, ElementRegionBase & )
  {
    localIndex const modelIndex = constitutiveManager->GetSubGroups().getIndex( m_fluidModelNames[targetIndex] );
    regionFluidMap.emplace( er, modelIndex );
  } );

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex > > dofNumberAccessor =
    elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex > >( dofKey );

  FluxKernel::ElementView< arrayView1d< globalIndex const > > const & dofNumber = dofNumberAccessor.toViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dPres       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & gravCoef    = m_gravCoef;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dens        = m_density;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & mob         = m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dMob_dPres  = m_dMobility_dPres;

  ElementRegionManager::ConstitutiveRelationAccessor< SingleFluidBase > fluidModels =
    elemManager.ConstructFullConstitutiveAccessor< SingleFluidBase >( constitutiveManager );

  // use ArrayView to make capture by value easy in lambdas
  arrayView1d< real64 const > const & presFace      = faceManager.getReference< array1d< real64 > >( viewKeyStruct::boundaryFacePressureString );
  arrayView2d< real64 >       const & densFace      = faceManager.getReference< array2d< real64 > >( viewKeyStruct::boundaryFaceDensityString );
  arrayView2d< real64 >       const & viscFace      = faceManager.getReference< array2d< real64 > >( viewKeyStruct::boundaryFaceViscosityString );
  arrayView1d< real64 >       const & mobFace       = faceManager.getReference< array1d< real64 > >( viewKeyStruct::boundaryFaceMobilityString );
  arrayView1d< real64 const > const & gravCoefFace  = faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  dataRepository::Group const & sets = faceManager.sets();

  // first, evaluate BC to get primary field values (pressure)
//  fsManager->ApplyField(faceManager, viewKeyStruct::boundaryFacePressure, time + dt);
  fsManager.Apply( time_n + dt,
                   domain,
                   "faceManager",
                   viewKeyStruct::boundaryFacePressureString,
                   [&] ( FieldSpecificationBase const * const fs,
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * const targetGroup,
                         string const fieldName )
  {
    fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet, time_n + dt, targetGroup, fieldName );
  } );


  // call constitutive models to get dependent quantities needed for flux (density, viscosity)
  fsManager.Apply( time_n + dt,
                   domain,
                   "faceManager",
                   viewKeyStruct::boundaryFacePressureString,
                   [&] ( FieldSpecificationBase const * GEOSX_UNUSED_PARAM( bc ),
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * const,
                         string const & )
  {
    for( auto kf : targetSet )
    {
      // since we don't have models on faces yet, we take them from an adjacent cell
      integer ke;
      for( ke = 0; ke < 2; ++ke )
      {
        if( elemRegionList[kf][ke] >= 0 && regionFluidMap.count( elemRegionList[kf][ke] ) > 0 )
        {
          break;
        }
      }
      GEOSX_ERROR_IF( ke > 1, "Face not adjacent to target regions: " << kf );
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];

      real64 dummy; // don't need derivatives on faces

      SingleFluidBase * const fluid = fluidModels[er][esr][regionFluidMap[er]];
      fluid->Compute( presFace[kf], densFace[kf][0], dummy, viscFace[kf][0], dummy );
    }

    MobilityKernel::Launch( targetSet, densFace, viscFace, mobFace );
  } );

  // *** assembly loop ***

  constexpr localIndex numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  real64 densWeight[numElems] = { 0.5, 0.5 };

  fsManager.Apply( time_n + dt,
                   domain,
                   "faceManager",
                   viewKeyStruct::boundaryFacePressureString,
                   [&] ( FieldSpecificationBase const * GEOSX_UNUSED_PARAM( bc ),
                         string const & setName,
                         SortedArrayView< localIndex const > const &,
                         Group * const,
                         string const & )
  {
    if( !sets.hasWrapper( setName ) || !fluxApprox->hasBoundaryStencil( setName ))
      return;

    FluxApproximationBase::BoundaryStencil const & stencil = fluxApprox->getBoundaryStencil( setName );
    ArrayOfArraysView< FluxApproximationBase::BoundaryStencil::Entry const, true > const & connections = stencil.getConnections();

    forAll< serialPolicy >( connections.size(), [=] ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray( iconn );

      stackArray1d< globalIndex, maxStencilSize > dofColIndices( stencilSize );

      stackArray1d< real64, numElems > mobility( numElems );
      stackArray1d< real64, numElems > dMobility_dP( numElems );
      stackArray1d< real64, maxStencilSize > dDensMean_dP( stencilSize );
      stackArray1d< real64, maxStencilSize > dFlux_dP( stencilSize );
      stackArray1d< real64, maxStencilSize > localFluxJacobian( stencilSize );

      // clear working arrays
      dDensMean_dP = 0.0;

      // calculate quantities on primary connected points
      real64 densMean = 0.0;
      globalIndex eqnRowIndex = -1;
      localIndex cell_order = -1;

      for( localIndex i = 0; i < numElems; ++i )
      {
        PointDescriptor const & point = connections( iconn, i ).index;

        real64 density = 0, dDens_dP = 0;
        switch( point.tag )
        {
          case PointDescriptor::Tag::CELL:
            {
              localIndex const er  = point.cellIndex.region;
              localIndex const esr = point.cellIndex.subRegion;
              localIndex const ei  = point.cellIndex.index;

              eqnRowIndex = dofNumber[er][esr][ei];

              density  = dens[er][esr][ei][0];
              dDens_dP = dDens_dPres[er][esr][ei][0];

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
            GEOSX_ERROR( "Unsupported point type in stencil" );
        }

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      dofColIndices = -1;
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        FluxApproximationBase::BoundaryStencil::Entry const & entry = connections( iconn, i );
        PointDescriptor const & point = entry.index;

        real64 pressure = 0.0, gravD = 0.0;
        switch( point.tag )
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
            GEOSX_ERROR( "Unsupported point type in stencil" );
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
      for( localIndex ke = 0; ke < stencilSize; ++ke )
        dFlux_dP[ke] *= mobility[k_up];
      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      integer sign = (cell_order == 0 ? 1 : -1);
      real64 const localFlux =  dt * flux * sign;

      integer counter = 0;
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        // compress arrays, skipping face derivatives
        if( dofColIndices[ke] >= 0 )
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

namespace
{
typedef SinglePhaseFVM< SinglePhaseBase > NoProppant;
typedef SinglePhaseFVM< SinglePhaseProppantBase > Proppant;
REGISTER_CATALOG_ENTRY( SolverBase, NoProppant, std::string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, Proppant, std::string const &, Group * const )
}
} /* namespace geosx */
