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
 * @file CompositionalMultiphaseHybridFVM.cpp
 */

#include "CompositionalMultiphaseHybridFVM.hpp"

#include "mesh/DomainPartition.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ResidualNormKernel.hpp"
#include "mesh/CellElementSubRegion.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace isothermalCompositionalMultiphaseBaseKernels;
using namespace compositionalMultiphaseHybridFVMKernels;
using namespace mimeticInnerProduct;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_lengthTolerance( 0 )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM;
}

void CompositionalMultiphaseHybridFVM::registerDataOnMesh( Group & meshBodies )
{
  // 1) Register the elem-centered data
  CompositionalMultiphaseBase::registerDataOnMesh( meshBodies );

  // 2) Register the face data
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();

    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressure changes

    faceManager.registerField< flow::facePressure_n >( getName() );

    // Register the face data for global component fraction
    faceManager.registerField< flow::faceGlobalCompFraction >( getName() );

    // Register the face data for temperature
    faceManager.registerField< flow::faceTemperature >( getName() );



    // Register the bc face data for pressure
    faceManager.registerField< flow::bcPressure >( getName() );

    // Register the bc face data for global component fraction
    faceManager.registerField< flow::bcGlobalCompFraction >( getName() );

    // Register the bc face data for temperature
    faceManager.registerField< flow::bcTemperature >( getName() );

    // Register face-based constitutive properties for BC faces
    // These will store fluid properties evaluated at BC conditions
    faceManager.registerField< flow::facePhaseMobility >( getName() ).
      reference().resizeDimension< 1 >( m_numPhases );

    faceManager.registerField< flow::facePhaseMassDensity >( getName() ).
      reference().resizeDimension< 1 >( m_numPhases );

    faceManager.registerField< flow::facePhaseCompFraction >( getName() ).
      reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

    // Register boundary face indicator (1 for boundary faces with Dirichlet BCs, 0 for interior)
    // Used to skip flux continuity constraints for boundary faces
    faceManager.registerField< flow::isBoundaryFace >( getName() );

    // auxiliary data for the buoyancy coefficient
    faceManager.registerField< flow::mimGravityCoefficient >( getName() );
  } );
}

void CompositionalMultiphaseHybridFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  GEOS_THROW_IF( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ),
                 getCatalogName() << " " << getDataContext() <<
                 ": the HybridMimeticDiscretization must be selected with CompositionalMultiphaseHybridFVM",
                 InputError, getDataContext() );

  GEOS_THROW_IF( m_hasCapPressure,
                 getCatalogName() << " " << getDataContext() <<
                 ": capillary pressure is not yet supported by CompositionalMultiphaseHybridFVM",
                 InputError, getDataContext() );
}

void CompositionalMultiphaseHybridFVM::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
  if( dynamicCast< QuasiRTInnerProduct const * >( &mimeticInnerProductBase )  ||
      dynamicCast< SimpleInnerProduct const * >( &mimeticInnerProductBase ) )
  {
    GEOS_ERROR( getCatalogName() << " " << getDataContext() <<
                "The QuasiRT, and Simple inner products are only available in SinglePhaseHybridFVM" );
  }

  m_lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * 1e-8;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager & faceManager = mesh.getFaceManager();

    CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups();

    // in the flux kernel, we need to make sure that we act only on the target regions
    // for that, we need the following region filter
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
    }

    // check that multipliers are stricly larger than 0, which would work with SinglePhaseFVM, but not with SinglePhaseHybridFVM.
    // To deal with a 0 multiplier, we would just have to skip the corresponding face in the FluxKernel
    arrayView1d< real64 const > const & transMultiplier = faceManager.getField< flow::transMultiplier >();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      minVal.min( transMultiplier[iface] );
    } );

    GEOS_THROW_IF( minVal.get() <= 0.0,
                   getCatalogName() << " " << getDataContext() <<
                   ": the transmissibility multipliers used in SinglePhaseHybridFVM must be strictly larger than 0.0",
                   std::runtime_error );

    // Initialize face-based constitutive property arrays to zero to prevent uninitialized memory usage on GPU
    arrayView2d< real64, compflow::USD_PHASE > facePhaseMob = faceManager.getField< flow::facePhaseMobility >();
    arrayView2d< real64, compflow::USD_PHASE > facePhaseMassDens = faceManager.getField< flow::facePhaseMassDensity >();
    arrayView3d< real64, compflow::USD_PHASE_COMP > facePhaseCompFrac = faceManager.getField< flow::facePhaseCompFraction >();

    localIndex const numFaces = faceManager.size();
    forAll< parallelDevicePolicy<> >( numFaces, [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      for( integer ip = 0; ip < facePhaseMob.size( 1 ); ++ip )
      {
        facePhaseMob[iface][ip] = 0.0;
        facePhaseMassDens[iface][ip] = 0.0;
        for( integer ic = 0; ic < facePhaseCompFrac.size( 2 ); ++ic )
        {
          facePhaseCompFrac[iface][ip][ic] = 0.0;
        }
      }
    } );

    // Mark boundary faces (faces with Dirichlet BCs) to skip flux continuity constraint
    // Initialize all faces as interior (0), then mark boundary faces (1)
    arrayView1d< integer > const isBoundaryFaceView = faceManager.getReference< array1d< integer > >( flow::isBoundaryFace::key() );
    // isBoundaryFaceView is default-initialized to zero

    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
    fsManager.forSubGroups< FieldSpecificationBase >( [&]( FieldSpecificationBase const & fs )
    {
      string const & fieldName = fs.getFieldName();
      if( fieldName == flow::bcPressure::key() ||
          fieldName == flow::bcGlobalCompFraction::key() ||
          fieldName == flow::bcTemperature::key() )
      {
        for( string const & setName : fs.getSetNames() )
        {
          SortedArrayView< localIndex const > const targetSet = faceManager.getSet( setName ).toViewConst();
          forAll< serialPolicy >( targetSet.size(), [=]( localIndex const i )
          {
            isBoundaryFaceView[targetSet[i]] = 1;
          } );
        }
      }
    } );

    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOS_LOG_RANK_0( getCatalogName() << " " << getDataContext() << ": An aquifer boundary condition named " <<
                       bc.getName() << " was requested in the XML file. \n" <<
                       "This type of boundary condition is not yet supported by CompositionalMultiphaseHybridFVM and will be ignored" );
    } );
  } );

}

void CompositionalMultiphaseHybridFVM::precomputeData( MeshLevel & mesh, string_array const & regionNames )
{
  FlowSolverBase::precomputeData( mesh, regionNames );

  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefNumerator;
  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefDenominator;
  mimFaceGravCoefNumerator.resize( faceManager.size() );
  mimFaceGravCoefDenominator.resize( faceManager.size() );

  // node data

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getField< flow::transMultiplier >();

  arrayView1d< real64 > const mimFaceGravCoef =
    faceManager.getField< flow::mimGravityCoefficient >();

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  real64 const lengthTolerance = m_lengthTolerance;

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & subRegion )
  {
    arrayView2d< real64 const > const & elemCenter =
      subRegion.template getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
    string const & permModelName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
    arrayView3d< real64 const > const & elemPerm =
      getConstitutiveModel< PermeabilityBase >( subRegion, permModelName ).permeability();
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.template getReference< array1d< real64 > >( flow::gravityCoefficient::key() );
    arrayView1d< real64 const > const & elemVolume = subRegion.getElementVolume();
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    // here we precompute some quantities (mimFaceFracCoef) used in the FluxKernel to assemble the one-sided gravity term in the transport
    // scheme
    // This one-sided gravity term is currently always treated with TPFA, as in MRST.
    // In the future, I will change that (here and in the FluxKernel) to have a consistent inner product for the gravity term as well
    compositionalMultiphaseHybridFVMKernels::
      simpleKernelLaunchSelector< PrecomputeKernel,
                                  mimeticInnerProduct::TPFAInnerProduct >( subRegion.numFacesPerElement(),
                                                                           subRegion.size(),
                                                                           faceManager.size(),
                                                                           nodePosition,
                                                                           faceToNodes,
                                                                           elemCenter,
                                                                           elemVolume,
                                                                           elemPerm,
                                                                           elemGravCoef,
                                                                           elemToFaces,
                                                                           transMultiplier,
                                                                           lengthTolerance,
                                                                           mimFaceGravCoefNumerator.toView(),
                                                                           mimFaceGravCoefDenominator.toView(),
                                                                           mimFaceGravCoef );

  } );

}

void CompositionalMultiphaseHybridFVM::implicitStepSetup( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // setup the elem-centered fields
  CompositionalMultiphaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & facePres_n =
      faceManager.getField< flow::facePressure_n >();
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< flow::facePressure >();
    facePres_n.setValues< parallelDevicePolicy<> >( facePres );
  } );

}


void CompositionalMultiphaseHybridFVM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                  DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // for the volume balance equation, disable global coupling
  // this equation is purely local (not coupled to neighbors or other physics)
  dofManager.disableGlobalCouplingForEquation( viewKeyStruct::elemDofFieldString(),
                                               m_numComponents );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::faceDofFieldString(),
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );

}


void CompositionalMultiphaseHybridFVM::assembleFluxTerms( real64 const dt,
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();


    // get the face-based DOF numbers for the assembly
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

    // Get boundary face indicator (initialized during initializePostInitialConditionsPreSubGroups)
    arrayView1d< integer const > const isBoundaryFaceView =
      faceManager.getField< flow::isBoundaryFace >();

    // get the element dof numbers for the assembly
    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    // get the face-centered pressures
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< flow::facePressure >();

    // get the face-centered depth
    arrayView1d< real64 const > const & faceGravCoef =
      faceManager.getField< flow::gravityCoefficient >();
    arrayView1d< real64 const > const & mimFaceGravCoef =
      faceManager.getField< flow::mimGravityCoefficient >();

    // get the face-centered transMultiplier
    arrayView1d< real64 const > const & transMultiplier =
      faceManager.getField< flow::transMultiplier >();

    // get the face-to-nodes connectivity for the transmissibility calculation
    ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

    arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemList          = faceManager.elementList().toViewConst();


    // tolerance for transmissibility calculation
    real64 const lengthTolerance = m_lengthTolerance;

    FluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    FluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      PermeabilityBase const & permeabilityModel =
        getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() ) );

      mimeticInnerProductReducedDispatch( mimeticInnerProductBase,
                                          [&] ( auto const mimeticInnerProduct )
      {
        using IP_TYPE = std::remove_const_t< TYPEOFREF( mimeticInnerProduct ) >;
        kernelLaunchSelector< FluxKernel,
                              IP_TYPE >( subRegion.numFacesPerElement(),
                                         m_numComponents, m_numPhases,
                                         er, esr, subRegion,
                                         permeabilityModel,
                                         m_regionFilter.toViewConst(),
                                         nodePosition,
                                         elemRegionList,
                                         elemSubRegionList,
                                         elemList,
                                         faceToNodes,
                                         faceDofNumber,
                                         faceGhostRank,
                                         isBoundaryFaceView,
                                         facePres,
                                         faceGravCoef,
                                         mimFaceGravCoef,
                                         transMultiplier,
                                         compFlowAccessors.get( flow::phaseMobility{} ),
                                         compFlowAccessors.get( flow::dPhaseMobility{} ),
                                         compFlowAccessors.get( flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseMassDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} ),
                                         elemDofNumber.toNestedViewConst(),
                                         dofManager.rankOffset(),
                                         lengthTolerance,
                                         dt,
                                         m_useTotalMassEquation,
                                         localMatrix,
                                         localRhs );

      } );
    } );

  } );
}

void CompositionalMultiphaseHybridFVM::assembleStabilizedFluxTerms( real64 const dt,
                                                                    DomainPartition const & domain,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs ) const
{
  // stab not implemented
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "Stabilized flux not available for this flow solver" );
}

real64 CompositionalMultiphaseHybridFVM::scalingForSystemSolution( DomainPartition & domain,
                                                                   DofManager const & dofManager,
                                                                   arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  real64 scalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< flow::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< flow::pressureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< flow::globalCompDensityScalingFactor >();
      auto const subRegionData =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionScalingKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution );

      scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
    } );

    FaceManager const & faceManager = mesh.getFaceManager();

    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    real64 const maxRelativePresChange = m_maxRelativePresChange;

    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< real64 const > const & facePressure =
      faceManager.getField< flow::facePressure >();
    globalIndex const rankOffset = dofManager.rankOffset();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minFaceVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
      {
        real64 const facePres = facePressure[iface];
        real64 const absPresChange = LvArray::math::abs( localSolution[faceDofNumber[iface] - rankOffset] );
        if( facePres > isothermalCompositionalMultiphaseBaseKernels::minDensForDivision )
        {
          real64 const relativePresChange = absPresChange / facePres;
          if( relativePresChange > maxRelativePresChange )
          {
            minFaceVal.min( maxRelativePresChange / relativePresChange );
          }
        }
      }
    } );

    if( minFaceVal.get() < scalingFactor )
    {
      scalingFactor = minFaceVal.get();
    }
  } );
  return LvArray::math::max( MpiWrapper::min( scalingFactor ), m_minScalingFactor );
}

bool CompositionalMultiphaseHybridFVM::checkSystemSolution( DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer localCheck = 1;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    // check cell-centered variables for each region
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure =
        subRegion.getField< flow::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getField< flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< flow::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< flow::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< flow::globalCompDensityScalingFactor >();
      // check that pressure and component densities are non-negative
      auto const subRegionData =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                     m_allowNegativePressure,
                                                     compositionalMultiphaseUtilities::ScalingType::Global,
                                                     scalingFactor,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     elemDofKey,
                                                     subRegion,
                                                     localSolution );

      localCheck = std::min( localCheck, subRegionData.localMinVal );
    } );
  } );

  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseHybridFVM::applyBoundaryConditions( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Apply face-based Dirichlet boundary conditions
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

void CompositionalMultiphaseHybridFVM::applyFaceDirichletBC( real64 const time_n,
                                                             real64 const dt,
                                                             DofManager const & dofManager,
                                                             DomainPartition & domain,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  using namespace compositionalMultiphaseHybridFVMKernels;
  using namespace mimeticInnerProduct;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  // Get the inner product type from the discretization
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  string const & innerProductType = hmDiscretization.getReference< string >( HybridMimeticDiscretization::viewKeyStruct::innerProductTypeString() );

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Log message for Dirichlet BC application
  static char const faceBcLogMessage[] =
    "CompositionalMultiphaseHybridFVM {}: at time {}s, "
    "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
    "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
    "\nThe total number of target faces (including ghost faces) is {}."
    "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";

  // Apply Dirichlet BC by computing boundary fluxes
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    globalIndex const rankOffset = dofManager.rankOffset();

    // Get node positions
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePosition = nodeManager.referencePosition();

    // Get face data
    ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();

    arrayView1d< real64 const > const transMultiplier = faceManager.getReference< array1d< real64 > >( fields::flow::transMultiplier::key() );

    // Apply boundary values to face fields first
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::bcPressure::key(), flow::facePressure::key() );

    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::bcGlobalCompFraction::key(), flow::faceGlobalCompFraction::key() );

    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    flow::bcTemperature::key(), flow::faceTemperature::key() );

    // Get face boundary values
    arrayView1d< real64 const > const facePres = faceManager.getField< fields::flow::facePressure >();
    arrayView1d< real64 const > const faceTemp = faceManager.getField< fields::flow::faceTemperature >();
    arrayView2d< real64 const, compflow::USD_COMP > const faceCompFrac = faceManager.getField< fields::flow::faceGlobalCompFraction >();
    arrayView1d< real64 const > const faceGravCoef = faceManager.getField< fields::flow::gravityCoefficient >();

    // Loop over regions and apply Dirichlet flux kernel
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const erIndex,
                                                                                CellElementSubRegion & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeabilityModel = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase & relperm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relPermName );

      real64 const lengthTolerance = m_lengthTolerance;

      // Get all face sets that have Dirichlet BCs
      std::set< string > bcFaceSets;
      fsManager.forSubGroups< FieldSpecificationBase >( [&]( FieldSpecificationBase const & fs )
      {
        string const & fieldName = fs.getFieldName();
        if( fieldName == flow::bcPressure::key() ||
            fieldName == flow::bcGlobalCompFraction::key() ||
            fieldName == flow::bcTemperature::key() )
        {
          for( string const & setName : fs.getSetNames() )
          {
            bcFaceSets.insert( setName );
          }
        }
      } );

      // Get writable face arrays for storing BC face properties
      arrayView2d< real64, compflow::USD_PHASE > facePhaseMob =
        faceManager.getField< flow::facePhaseMobility >();
      arrayView2d< real64, compflow::USD_PHASE > facePhaseMassDens =
        faceManager.getField< flow::facePhaseMassDensity >();
      arrayView3d< real64, compflow::USD_PHASE_COMP > facePhaseCompFrac =
        faceManager.getField< flow::facePhaseCompFraction >();

      // Move arrays to host memory before evaluateBCFaceProperties runs with serialPolicy
      facePhaseMob.move( hostMemorySpace, true );
      facePhaseMassDens.move( hostMemorySpace, true );
      facePhaseCompFrac.move( hostMemorySpace, true );

      // Evaluate constitutive properties at BC face conditions for each face set
      for( string const & setName : bcFaceSets )
      {
        SortedArrayView< localIndex const > const targetSet = faceManager.getSet( setName ).toViewConst();
        if( targetSet.size() == 0 )
          continue;

        // Call evaluateBCFaceProperties to compute face properties at BC conditions
        constitutive::constitutiveComponentUpdatePassThru( fluid, m_numComponents, [&]( auto & fluidWrapper, auto NC )
        {
          GEOS_UNUSED_VAR( fluidWrapper );
          integer constexpr NUM_COMP = NC();

          auto evaluateWithPhases = [&]( auto NP_VALUE )
          {
            integer constexpr NUM_PHASES = decltype( NP_VALUE )::value;

            compositionalMultiphaseHybridFVMKernels::evaluateBCFaceProperties< NUM_COMP, NUM_PHASES >
              ( m_numPhases,
              targetSet,
              facePres,
              faceTemp,
              faceCompFrac,
              elemRegionList,
              elemSubRegionList,
              elemList,
              erIndex,
              subRegion.getIndexInParent(),
              fluid,
              relperm,
              facePhaseMob,
              facePhaseMassDens,
              facePhaseCompFrac );
          };

          if( m_numPhases == 2 )
          {
            evaluateWithPhases( std::integral_constant< integer, 2 >() );
          }
          else if( m_numPhases == 3 )
          {
            evaluateWithPhases( std::integral_constant< integer, 3 >() );
          }
        } );
      }

      // CRITICAL: Move the SAME array views that were just modified back to device memory
      // Don't get new views - use the views that evaluateBCFaceProperties actually modified
      // evaluateBCFaceProperties uses serialPolicy (host), DirichletFluxKernel uses parallelDevicePolicy (device)
      facePhaseMob.move( parallelDeviceMemorySpace, false );
      facePhaseMassDens.move( parallelDeviceMemorySpace, false );
      facePhaseCompFrac.move( parallelDeviceMemorySpace, false );

      // Get const views to the face properties for use in DirichletFluxKernel
      arrayView2d< real64 const, compflow::USD_PHASE > const facePhaseMobField =
        faceManager.getField< flow::facePhaseMobility >();
      arrayView2d< real64 const, compflow::USD_PHASE > const facePhaseMassDensField =
        faceManager.getField< flow::facePhaseMassDensity >();
      arrayView3d< real64 const, compflow::USD_PHASE_COMP > const facePhaseCompFracField =
        faceManager.getField< flow::facePhaseCompFraction >();

      // Apply Dirichlet boundary fluxes for each face set using DirichletFluxKernel
      for( string const & setName : bcFaceSets )
      {
        SortedArrayView< localIndex const > const targetSet = faceManager.getSet( setName ).toViewConst();
        if( targetSet.size() == 0 )
          continue;

        // Launch the Dirichlet flux kernel with compile-time dispatch
        constitutive::constitutiveComponentUpdatePassThru( fluid, m_numComponents, [&]( auto & fluidWrapper, auto NC )
        {
          GEOS_UNUSED_VAR( fluidWrapper );
          integer constexpr NUM_COMP = NC();

          typename DirichletFluxKernel::CompFlowAccessors compFlowAccessors( elemManager, getName() );
          typename DirichletFluxKernel::MultiFluidAccessors multiFluidAccessors( elemManager, getName() );
          ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumberAccessor =
            elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );

          localIndex const numFacesPerElement = subRegion.numFacesPerElement();

          auto launchWithPhases = [&]( auto NP_VALUE )
          {
            integer constexpr NP = decltype( NP_VALUE )::value;

            auto launchKernel = [&]( auto IP_TYPE_WRAPPER, auto NF_VALUE )
            {
              using IP_TYPE = decltype( IP_TYPE_WRAPPER );
              integer constexpr NF = decltype( NF_VALUE )::value;

              DirichletFluxKernel::launch< NF, NUM_COMP, NP, IP_TYPE >
                ( m_numPhases,
                erIndex,
                subRegion.getIndexInParent(),
                subRegion,
                permeabilityModel,
                targetSet,
                facePres,
                faceTemp,
                facePhaseMobField,
                facePhaseMassDensField,
                facePhaseCompFracField,
                nodePosition,
                faceToNodes,
                elemRegionList,
                elemSubRegionList,
                elemList,
                faceDofNumber,
                faceGravCoef,
                transMultiplier,
                lengthTolerance,
                dt,
                rankOffset,
                m_useTotalMassEquation,
                compFlowAccessors,
                multiFluidAccessors,
                elemDofNumberAccessor.toNestedViewConst(),
                localMatrix,
                localRhs );
            };

            // Helper lambda to dispatch on number of faces per element
            auto launchForFaces = [&]( auto IP_TAG )
            {
              if( numFacesPerElement == 4 )
                launchKernel( IP_TAG, std::integral_constant< integer, 4 >{} );
              else if( numFacesPerElement == 5 )
                launchKernel( IP_TAG, std::integral_constant< integer, 5 >{} );
              else if( numFacesPerElement == 6 )
                launchKernel( IP_TAG, std::integral_constant< integer, 6 >{} );
              else if( numFacesPerElement == 7 )
                launchKernel( IP_TAG, std::integral_constant< integer, 7 >{} );
              else if( numFacesPerElement == 8 )
                launchKernel( IP_TAG, std::integral_constant< integer, 8 >{} );
              else if( numFacesPerElement == 9 )
                launchKernel( IP_TAG, std::integral_constant< integer, 9 >{} );
              else if( numFacesPerElement == 10 )
                launchKernel( IP_TAG, std::integral_constant< integer, 10 >{} );
              else if( numFacesPerElement == 11 )
                launchKernel( IP_TAG, std::integral_constant< integer, 11 >{} );
              else if( numFacesPerElement == 12 )
                launchKernel( IP_TAG, std::integral_constant< integer, 12 >{} );
              else if( numFacesPerElement == 13 )
                launchKernel( IP_TAG, std::integral_constant< integer, 13 >{} );
              else
                GEOS_ERROR( "Unsupported number of faces per element: " << numFacesPerElement );
            };

            // Inner-product selection
            if( innerProductType == "TPFA" )
            {
              launchForFaces( TPFAInnerProduct{} );
            }
            else if( innerProductType == "quasiTPFA" )
            {
              launchForFaces( QuasiTPFAInnerProduct{} );
            }
            else if( innerProductType == "beiraoDaVeigaLipnikovManzini" )
            {
              launchForFaces( BdVLMInnerProduct{} );
            }
            else
            {
              GEOS_ERROR( "Unsupported inner product type: " << innerProductType );
            }
          };

          if( m_numPhases == 2 )
            launchWithPhases( std::integral_constant< integer, 2 >{} );
          else if( m_numPhases == 3 )
            launchWithPhases( std::integral_constant< integer, 3 >{} );
          else
            GEOS_ERROR( "Unsupported number of phases: " << m_numPhases );
        } );
      }
    } );

    // Keep strong enforcement that the face pressure equals the informed bcPressure value (original behavior)
    arrayView1d< real64 const > const presFace =
      faceManager.getField< flow::facePressure >();
    arrayView1d< real64 const > const presFaceBC =
      faceManager.getField< flow::bcPressure >();

    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    flow::bcPressure::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager & targetGroup,
                                          string const & )
    {
      GEOS_UNUSED_VAR( setName );
      // Using the field specification functions to apply the boundary conditions to the system
      fs.applyFieldValue< FieldSpecificationEqual,
                          parallelDevicePolicy<> >( targetSet,
                                                    time_n + dt,
                                                    targetGroup,
                                                    flow::bcPressure::key() );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {

        localIndex const kf = targetSet[a];
        if( faceGhostRank[kf] >= 0 )
        {
          return;
        }

        // Get the dof number of this face
        globalIndex const dofIndex = faceDofNumber[kf];
        localIndex const localRow = dofIndex - rankOffset;

        // ENFORCE DIRICHLET CONSTRAINT: x_i = x_spec
        // Mathematical procedure to enforce prescribed value in Ax = b:
        // 1. For row i: Zero all entries except diagonal, set A[i,i] = 1, set b[i] = x_spec
        // 2. For all other rows k: subtract A[k,i] * x_spec from b[k], then set A[k,i] = 0
        //
        // Note: Step 2 (removing column influence) should ideally be done, but in practice
        // for boundary face DOFs in hybrid FVM, the column entries from interior faces to
        // boundary faces are typically zero or minimal because boundary faces don't strongly
        // couple to interior equations. The strong coupling is boundary->interior (via fluxes).
        // For exact enforcement in non-trivial cases, column zeroing would be needed.

        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( localRow );
          arraySlice1d< real64 > const entries = localMatrix.getEntries( localRow );
          localIndex const numEntries = localMatrix.numNonZeros( localRow );

          // Step 1: Zero out all row entries and set diagonal to 1
          for( localIndex j = 0; j < numEntries; ++j )
          {
            if( columns[j] == dofIndex )
            {
              entries[j] = 1.0;  // Set diagonal to 1
            }
            else
            {
              entries[j] = 0.0;  // Zero out off-diagonal entries
            }
          }

          // Set RHS to the prescribed boundary value (absolute value)
          localRhs[localRow] = presFace[kf] - presFaceBC[kf];
        }

      } );
    } );

  } );
}

void CompositionalMultiphaseHybridFVM::applyAquiferBC( real64 const time,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );

}

void CompositionalMultiphaseHybridFVM::saveAquiferConvergedState( real64 const & time,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt, domain );
}


real64 CompositionalMultiphaseHybridFVM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                                real64 const & dt,
                                                                DomainPartition const & domain,
                                                                DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 localResidualNorm = 0.0;
  real64 localResidualNormalizer = 0.0;

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  // local residual
  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // here we compute the cell-centered residual norm in the derived class
    // to avoid duplicating a synchronization point

    // get a view into local residual vector

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[2]{};
      real64 subRegionResidualNormalizer[2]{};

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1.1: compute the norm in the subRegion

      isothermalCompositionalMultiphaseBaseKernels::
        ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( normType,
                                                   numFluidComponents(),
                                                   rankOffset,
                                                   elemDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   solid,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm,
                                                   subRegionResidualNormalizer );

      // step 1.2: reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        // take max between mass and volume residual
        subRegionResidualNorm[0] = LvArray::math::max( subRegionResidualNorm[0], subRegionResidualNorm[1] );
        if( subRegionResidualNorm[0] > localResidualNorm )
        {
          localResidualNorm = subRegionResidualNorm[0];
        }
      }
      else
      {
        // sum up mass and volume residual
        subRegionResidualNorm[0] = subRegionResidualNorm[0] + subRegionResidualNorm[1];
        localResidualNorm += subRegionResidualNorm[0];
        localResidualNormalizer += subRegionResidualNormalizer[0];
      }
    } );

    // step 2: compute the residual for the face-based constraints

    real64 faceResidualNorm[1]{};
    real64 faceResidualNormalizer[1]{};

    // step 2.1: compute the norm for the local faces

    compositionalMultiphaseHybridFVMKernels::
      ResidualNormKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( normType,
                                                 rankOffset,
                                                 faceDofKey,
                                                 localRhs,
                                                 m_regionFilter.toViewConst(),
                                                 getName(),
                                                 elemManager,
                                                 faceManager,
                                                 dt,
                                                 m_nonlinearSolverParameters.m_minNormalizer,
                                                 faceResidualNorm,
                                                 faceResidualNormalizer );

    // step 2.2: reduction across meshBodies/regions/subRegions

    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      if( faceResidualNorm[0] > localResidualNorm )
      {
        localResidualNorm = faceResidualNorm[0];
      }
    }
    else
    {
      localResidualNorm += faceResidualNorm[0];
      localResidualNormalizer += faceResidualNormalizer[0];
    }
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm, residualNorm );
  }
  else
  {
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm, localResidualNormalizer, residualNorm );
  }

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm ));

  getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), residualNorm );

  return residualNorm;
}

void CompositionalMultiphaseHybridFVM::applySystemSolution( DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor,
                                                            real64 const GEOS_UNUSED_PARAM( dt ),
                                                            DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // 1. apply the elem-based update
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               flow::globalCompDensity::key(),
                               scalingFactor,
                               ~pressureMask );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::faceDofFieldString(),
                               flow::facePressure::key(),
                               scalingFactor );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    {
      fieldsToBeSync.addElementFields( { flow::pressure::key(),
                                         flow::globalCompDensity::key() },
                                       regionNames );

      fieldsToBeSync.addFields( FieldLocation::Face, { flow::facePressure::key() } );
    };

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}


void CompositionalMultiphaseHybridFVM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // 1. Reset the cell-centered fields
  CompositionalMultiphaseBase::resetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & facePres_n =
      faceManager.getField< flow::facePressure_n >();
    arrayView1d< real64 > const & facePres =
      faceManager.getField< flow::facePressure >();
    facePres.setValues< parallelDevicePolicy<> >( facePres_n );
  } );
}

void CompositionalMultiphaseHybridFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  MultiFluidBase const & fluid =
    getConstitutiveModel< MultiFluidBase >( dataGroup,
                                            dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  RelativePermeabilityBase const & relperm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup,
                                                      dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() ) );

  compositionalMultiphaseHybridFVMKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geos */
