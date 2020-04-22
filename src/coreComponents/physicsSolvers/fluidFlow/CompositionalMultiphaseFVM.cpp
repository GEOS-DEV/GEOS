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
 * @file CompositionalMultiphaseFVM.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "dataRepository/Group.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseFVMKernels;

CompositionalMultiphaseFVM::CompositionalMultiphaseFVM( const string & name,
                                                        Group * const parent )
  :
  CompositionalMultiphaseBase( name, parent )
{}


void CompositionalMultiphaseFVM::SetupDofs( DomainPartition const * const domain,
                                            DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString, fluxApprox );
}


void CompositionalMultiphaseFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                    real64 const dt,
                                                    DomainPartition const * const domain,
                                                    DofManager const * const dofManager,
                                                    ParallelMatrix * const matrix,
                                                    ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d< globalIndex >,
                                                          arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementView< arrayView1d< globalIndex const > > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementView< arrayView1d< real64 const > > const & pres                = m_pressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & dPres               = m_deltaPressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & gravCoef            = m_gravCoef.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & phaseMob            = m_phaseMob.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres     = m_dPhaseMob_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp     = m_dPhaseMob_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp = m_dPhaseVolFrac_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens = m_dCompFrac_dCompDens.toViewConst();

  FluxKernel::ElementView< arrayView3d< real64 const > > const & phaseDens                   = m_phaseDens.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres            = m_dPhaseDens_dPres.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp            = m_dPhaseDens_dComp.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & phaseCompFrac               = m_phaseCompFrac.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres        = m_dPhaseCompFrac_dPres.toViewConst();
  FluxKernel::ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp        = m_dPhaseCompFrac_dComp.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & phaseCapPres                = m_phaseCapPressure.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseCapPres_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac.toViewConst();

  localIndex constexpr numElems   = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex constexpr maxSize1 = numElems * maxNumComp;
  localIndex constexpr maxSize2 = maxStencil * maxNumDof;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  integer const capPressureFlag = m_capPressureFlag;

  fluxApprox->forAllStencils( [&] ( auto const & stencil )
  {
    typedef TYPEOFREF( stencil ) STENCIL_TYPE;
    typename STENCIL_TYPE::IndexContainerViewConstType const & eri = stencil.getElementRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & esri = stencil.getElementSubRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & ei = stencil.getElementIndices();
    typename STENCIL_TYPE::WeightContainerViewConstType const & weights = stencil.getWeights();

    forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
    {
      localIndex const stencilSize = stencil.stencilSize( iconn );

      // create local work arrays
      stackArray1d< globalIndex, maxSize1 > eqnRowIndices( numElems * NC );
      stackArray1d< globalIndex, maxSize2 > dofColIndices( stencilSize * NDOF );

      stackArray1d< real64, maxSize1 >            localFlux( numElems * NC );
      stackArray2d< real64, maxSize1 * maxSize2 > localFluxJacobian( numElems * NC, stencilSize * NDOF );

      FluxKernel::Compute( NC, NP,
                           stencilSize,
                           eri[iconn],
                           esri[iconn],
                           ei[iconn],
                           weights[iconn],
                           pres,
                           dPres,
                           gravCoef,
                           phaseMob,
                           dPhaseMob_dPres,
                           dPhaseMob_dComp,
                           dPhaseVolFrac_dPres,
                           dPhaseVolFrac_dComp,
                           dCompFrac_dCompDens,
                           phaseDens,
                           dPhaseDens_dPres,
                           dPhaseDens_dComp,
                           phaseCompFrac,
                           dPhaseCompFrac_dPres,
                           dPhaseCompFrac_dComp,
                           phaseCapPres,
                           dPhaseCapPres_dPhaseVolFrac,
                           capPressureFlag,
                           dt,
                           localFlux,
                           localFluxJacobian );

      // set equation indices for both connected cells
      for( localIndex i = 0; i < numElems; ++i )
      {
        globalIndex const offset = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          eqnRowIndices[i * NC + ic] = offset + ic;
        }
      }

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        globalIndex const offset = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];

        for( localIndex jdof = 0; jdof < NDOF; ++jdof )
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }
      }

      // TODO: apply equation/variable change transformation(s)

      // Add to global residual/jacobian
      rhs->add( eqnRowIndices.data(),
                localFlux.data(),
                numElems * NC );

      matrix->add( eqnRowIndices.data(),
                   dofColIndices.data(),
                   localFluxJacobian.data(),
                   numElems * NC,
                   stencilSize * NDOF );

    } );
  } );
}


void CompositionalMultiphaseFVM::ApplyBoundaryConditions( real64 const time_n,
                                                          real64 const dt,
                                                          DomainPartition * const domain,
                                                          DofManager const & dofManager,
                                                          ParallelMatrix & matrix,
                                                          ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit( time_n, dt, &dofManager, domain, &matrix, &rhs );

  // apply flux boundary conditions
  ApplySourceFluxBC( time_n, dt, &dofManager, domain, &matrix, &rhs );


  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void
CompositionalMultiphaseFVM::ApplySourceFluxBC( real64 const time,
                                               real64 const dt,
                                               DofManager const * const dofManager,
                                               DomainPartition * const domain,
                                               ParallelMatrix * const matrix,
                                               ParallelVector * const rhs )
{


  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  fsManager.Apply( time + dt, domain, "ElementRegions", FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
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

    fs->ApplyBoundaryConditionToSystem< FieldSpecificationAdd, LAInterface >( localSet,
                                                                              time + dt,
                                                                              dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              integer_conversion< int >( m_numDofPerCell ),
                                                                              *matrix,
                                                                              *rhs,
                                                                              [&] ( localIndex const GEOSX_UNUSED_PARAM( a )) -> real64
    {
      return 0;
    } );

  } );
}


void
CompositionalMultiphaseFVM::ApplyDirichletBC_implicit( real64 const time,
                                                       real64 const dt,
                                                       DofManager const * const dofManager,
                                                       DomainPartition * const domain,
                                                       ParallelMatrix * const matrix,
                                                       ParallelVector * const rhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const & setName,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * subRegion,
                        string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    string const & subRegionName = subRegion->getName();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting pressure boundary conditions on set " << setName );
    bcStatusMap[subRegionName][setName].resize( m_numComponents );
    bcStatusMap[subRegionName][setName] = false;

    // 1.1. Apply BC to set the field values
    fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                    time + dt,
                                                    subRegion,
                                                    viewKeyStruct::bcPressureString );
  } );

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::globalCompFractionString,
                   [&] ( FieldSpecificationBase const * const fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * subRegion,
                         string const & )
  {
    // 2.0. Check pressure and record composition bc application
    string const & subRegionName = subRegion->getName();
    localIndex const comp = fs->GetComponent();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[subRegionName][setName][comp] = true;

    // 2.1. Apply BC to set the field values
    fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                    time + dt,
                                                    subRegion,
                                                    viewKeyStruct::globalCompFractionString );
  } );

  // 2.3 Check consistency between composition BC applied to sets
  bool bcConsistent = true;
  for( auto const & bcStatusEntryOuter : bcStatusMap )
  {
    for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
    {
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        bcConsistent &= bcStatusEntryInner.second[ic];
        GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic <<
                          " on region " << bcStatusEntryOuter.first << ", set " << bcStatusEntryInner.first );
      }
    }
  }
  GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString,
                   [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                         string const & GEOSX_UNUSED_PARAM( setName ),
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * subRegion,
                         string const & )
  {
    // TODO: hack! Find a better way to get the fluid
    Group const * const region = subRegion->getParent()->getParent();
    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region->getName() ) ];
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( *subRegion, fluidName );

    arrayView1d< globalIndex const > const & dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const & pres      = subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dPres     = subRegion->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & bcPres    = subRegion->getReference< array1d< real64 > >( viewKeyStruct::bcPressureString );
    arrayView2d< real64 const > const & compFrac  = subRegion->getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
    arrayView2d< real64 const > const & compDens  = subRegion->getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dCompDens = subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    array1d< real64 > rhsContribution( targetSet.size() * m_numDofPerCell );
    array1d< globalIndex > dof( targetSet.size() * m_numDofPerCell );

    integer counter = 0;

    for( localIndex a : targetSet )
    {
      fluid.PointUpdate( bcPres[a], m_temperature, compFrac[a], a, 0 );

      dof[counter] = dofNumber[a];

      // 4.1. Apply pressure to the matrix
      FieldSpecificationEqual::SpecifyFieldValue< LAInterface >( dof[counter],
                                                                 *matrix,
                                                                 rhsContribution[counter],
                                                                 bcPres[a],
                                                                 pres[a] + dPres[a] );

      ++counter;

      // 4.2. For each component, apply target global density value
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        dof[counter] = dofNumber[a] + ic + 1;
        real64 const targetCompDens = totalDens[a][0] * compFrac[a][ic];

        FieldSpecificationEqual::SpecifyFieldValue< LAInterface >( dof[counter],
                                                                   *matrix,
                                                                   rhsContribution[counter],
                                                                   targetCompDens,
                                                                   compDens[a][ic] + dCompDens[a][ic] );

        ++counter;
      }

      // 4.3. Apply accumulated rhs values
      FieldSpecificationEqual::PrescribeRhsValues< LAInterface >( *rhs,
                                                                  counter,
                                                                  dof.data(),
                                                                  rhsContribution.data() );
    }
  } );
}

real64
CompositionalMultiphaseFVM::CalculateResidualNorm( DomainPartition const * const domain,
                                                   DofManager const & dofManager,
                                                   ParallelVector const & rhs )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  real64 const * localResidual = rhs.extractLocalVector();
  real64 localResidualNorm = 0.0;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase const & subRegion )
  {
    MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & refPoro =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        globalIndex const offset = dofNumber[ei];
        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          localIndex const lid = rhs.getLocalRowID( offset + idof );
          real64 const val = localResidual[lid] / ( totalDens[ei][0] * refPoro[ei] * volume[ei] );
          localResidualNorm += val * val;
        }
      }
    }
  } );

  // compute global residual norm
  realT globalResidualNorm;
  MpiWrapper::allReduce( &localResidualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX );

  return sqrt( globalResidualNorm ); 
}

bool
CompositionalMultiphaseFVM::CheckSystemSolution( DomainPartition const * const domain,
                                                 DofManager const & dofManager,
                                                 ParallelVector const & solution,
                                                 real64 const scalingFactor )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  real64 const * localSolution = solution.extractLocalVector();
  int localCheck = 1;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & pres = m_pressure[er][esr];
    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];
    arrayView2d< real64 const > const & compDens = m_globalCompDensity[er][esr];
    arrayView2d< real64 const > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        globalIndex const offset = dofNumber[ei];
        // extract solution and apply to dP
        {
          localIndex const lid = solution.getLocalRowID( offset );
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

          if( newPres < 0.0 )
          {
            localCheck = 0;
          }
        }

        for( localIndex ic = 0; ic < m_numComponents; ++ic )
        {
          localIndex const lid = solution.getLocalRowID( offset + ic + 1 );
          real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[lid];
          if( newDens < 0.0 )
          {
            localCheck = 0;
          }
        }
      }
    } );
  } );

  return MpiWrapper::Min( localCheck, MPI_COMM_GEOSX );
}

void
CompositionalMultiphaseFVM::ApplySystemSolution( DofManager const & dofManager,
                                                 ParallelVector const & solution,
                                                 real64 const scalingFactor,
                                                 DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerCell );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain->getNeighbors() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFVM, string const &, Group * const )
}// namespace geosx
