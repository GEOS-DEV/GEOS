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
 * @file SolidMechanicsMPM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_

#include "physicsSolvers/PhysicsSolverBase.hpp"

#include "common/TimingMacros.hpp"

#include "kernels/ExplicitMPM.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"

// #include "kernels/SolidMechanicsLagrangianFEMKernels.hpp"
// #include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

#include "events/mpmEvents/MPMEventManager.hpp"
#include "mesh/CohesiveZoneManager.hpp"

#include "physicsSolvers/solidMechanics/MPMSolverFields.hpp"
#include "physicsSolvers/solidMechanics/MPMEnums.hpp"

namespace geos
{

class SpatialPartition;


/**
 * @class SolidMechanicsMPM
 *
 * This class implements a material point method solution to the equations of motion.
 */
class SolidMechanicsMPM : public PhysicsSolverBase
{
public:

  /**
   * @enum TimeIntegrationOption
   *
   * The options for time integration
   */
  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,      //!< QuasiStatic
    ImplicitDynamic,  //!< ImplicitDynamic
    ExplicitDynamic   //!< ExplicitDynamic
  };

  /**
   * @enum UpdateMethodOption
   *
   * The options for time integration
   */
  enum class UpdateMethodOption : integer
  {
    FLIP,      //!< FLIP
    PIC,       //!< PIC
    XPIC,      //!< XPIC
    FMPM       //!< FMPM
  };

  /**
   * @enum BoundaryConditionOption
   *
   * The options for essential boundary conditions
   */
  enum struct BoundaryConditionOption : integer
  {
    Outflow,    //!<Outflow
    Symmetry,   //!<Symmetry
    Moving,     //!<Moving
    Contact     //!<Contact
  };

  /**
   * @enum SurfaceFlag
   *
   * The flags associated with different surface types
   */
  enum struct SurfaceFlag : integer
  {
    Interior,
    FullyDamaged,
    Surface,
    Cohesive,
    DamagedCohesive
  };

  /**
   * @enum ContactNormalTypeOption
   *
   * The options for contact gap correction
   */
  enum struct ContactNormalTypeOption : integer
  {
    Difference,
    MassWeighted,
    LargerMass,
    Mixed,
    Aligned,
    LogisticRegression
  };

  /**
   * @enum ContactGapCorrectionOption
   *
   * The options for contact gap correction
   */
  enum struct ContactGapCorrectionOption : integer
  {
    Simple,
    Implicit,
    Softened
  };

  /**
   * @enum AreaIntegrationOption
   *
   * The options for nodal area integration
   */
  enum struct AreaIntegrationOption : integer
  {
    BruteForce,
    Mesh
  };

  /**
   * @enum OverlapCorrectionOption
   *
   * The options for overlap correction
   */
  enum struct OverlapCorrectionOption : integer
  {
    Off,
    NormalForce,
    SPH,
  };


  /**
   * @enum CohesiveLawOption
   *
   * The options for cohesive laws
   */
  enum struct CohesiveLawOption : integer
  {
    Uncoupled,
    NeedlemanXu,
    Polymer
  };

  enum struct GPUSchemeOption : integer
  {
    Atomics,
    NoAtomics,
    RandomMix,
    MinimalAtomics,
    Reduction,
    Colors
  };

  enum struct NormalsAndPositionsMethodOption : integer
  {
    LogisticRegression,
    DFGAndVolumeIntegration
  };

  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsMPM( const string & name,
                     Group * const parent );


  SolidMechanicsMPM( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM( SolidMechanicsMPM && ) = default;

  SolidMechanicsMPM & operator=( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM & operator=( SolidMechanicsMPM && ) = delete;

  /**
   * destructor
   */
  virtual ~SolidMechanicsMPM() override {}

  /**
   * @return The string that may be used to generate a new instance from the PhysicsSolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanics_MPM"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final
  {
    // There should be nothing to update
    GEOS_UNUSED_VAR( domain );
  };

  /**@}*/

  /**
   * This method is called when its host event is triggered
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;


  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  void assemblyLaunch( DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs,
                       PARAMS && ... params );


  template< typename ... PARAMS >
  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 string_array const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );
                                 
  struct groupKeyStruct : PhysicsSolverBase::groupKeyStruct
  {
    static constexpr char const * mpmEventManagerString() { return "MPMEvents"; };
    static constexpr char const * cohesiveZoneManagerString() { return "CohesiveZoneManager"; };
  };

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }

    static constexpr char const * gridExternalForceString() { return "gridExternalForce"; }
    static constexpr char const * gridInternalForceString() { return "gridInternalForce"; }
    static constexpr char const * gridSurfaceTensionForceString() { return "gridSurfaceForceTension"; }
    static constexpr char const * gridDisplacementString() { return "gridDisplacement"; }
    static constexpr char const * gridCenterOfVolumeString() { return "gridCenterOfVolume"; }
    static constexpr char const * gridParticleMappedSurfaceNormalString() { return "gridParticleMappedSurfaceNormal"; }

    static constexpr char const * gridMassString() { return "gridMass"; }
    static constexpr char const * gridMaterialVolumeString() { return "gridMaterialVolume"; }
    static constexpr char const * gridVelocityString() { return "gridVelocity"; }
    static constexpr char const * gridDVelocityString() { return "gridDVelocity"; }
    static constexpr char const * gridMomentumString() { return "gridMomentum"; }
    static constexpr char const * gridAccelerationString() { return "gridAcceleration"; }
    static constexpr char const * gridContactForceString() { return "gridContactForce"; }
    static constexpr char const * gridDamageString() { return "gridDamage"; }
    static constexpr char const * gridDamageGradientString() { return "gridDamageGradient"; }
    static constexpr char const * gridMaxDamageString() { return "gridMaxDamage"; }

    static constexpr char const * gridSurfaceNormalWeightsString() { return "gridSurfaceNormalWeights"; }
    static constexpr char const * gridSurfaceNormalWeightNormalizationString() { return "gridSurfaceNormalWeightNormalization"; }
    static constexpr char const * gridSurfaceNormalString() { return "gridSurfaceNormal"; }
    static constexpr char const * gridSurfacePositionString() { return "gridSurfacePosition"; }
    static constexpr char const * gridSurfaceAreaString() { return "gridSurfaceArea"; }

    static constexpr char const * gridCenterOfMassString() { return "gridCenterOfMass"; }
    static constexpr char const * gridNormalStressString() { return "gridNormalStress"; }
    static constexpr char const * gridMassWeightedDamageString() { return "gridMassWeightedDamage"; }
    static constexpr char const * gridCohesiveNodeString() { return "gridCohesiveNode"; }
    static constexpr char const * gridReferenceAreaVectorString() { return "gridReferenceAreaVector"; }
    static constexpr char const * gridReferenceSurfacePositionString() { return "gridReferenceSurfacePosition"; }
    static constexpr char const * gridReferenceMaterialVolumeString() { return "gridReferenceMaterialVolume"; }

    static constexpr char const * gridBackgroundStressString() { return "gridBackgroundStress"; }

    static constexpr char const * gridSurfaceMassString() { return "gridSurfaceMass"; }
    static constexpr char const * gridSurfaceFieldMassString() { return "gridSurfaceFieldMass"; }
    static constexpr char const * gridExplicitSurfaceNormalString() { return "gridExplicitSurfaceNormal"; }
    static constexpr char const * gridMaxMappedParticleIDString() { return "gridMaxMappedParticleIDS"; }
    static constexpr char const * gridPrincipalExplicitSurfaceNormalString() { return "gridPrincipalExplicitSurfaceNormal"; }
    static constexpr char const * gridCohesiveFieldFlagString() { return "gridCohesiveFieldFlag"; }
    static constexpr char const * gridCohesiveAreaString() { return "gridCohesiveArea"; }
    static constexpr char const * gridCohesiveForceString() { return "gridCohesiveForce"; }

    static constexpr char const * gridHasBinderString() { return "gridHasBinder"; }

    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }

    static constexpr char const * gridVPlusString() { return "gridVPlus"; }
    static constexpr char const * gridDVPlusString() { return "gridDVPlus"; }

    static constexpr char const * gridBasedSurfaceNormalString() { return "gridBasedSurfaceNormal"; }
    static constexpr char const * gridBasedSurfacePositionString() { return "gridBasedSurfacePosition"; }

    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  void initialize( NodeManager & nodeManager,
                   ParticleManager & particleManager,
                   SpatialPartition & partition );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  localIndex partitionField( int numContactGroups,
                             int damageFieldPartitioning,
                             localIndex particleGroup,
                             arraySlice1d< real64 const > const particleDamageGradient,
                             arraySlice1d< real64 const > const particleSurfaceNormal,
                             arraySlice1d< real64 const > const gridDamageGradient );

  void triggerEvents( const real64 dt,
                      const real64 time_n,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  void checkEventCompletion( const real64 time_n );

  void performMaterialSwap( ParticleManager & particleManager,
                            string sourceRegionName,
                            string destinationRegionName );

  void resizeGrid( SpatialPartition & partition,
                   NodeManager & nodeManager,
                   real64 const dt );

  void syncGridFields( stdVector< std::string > const & fieldNames,
                       DomainPartition & domain,
                       NodeManager & nodeManager,
                       MeshLevel & mesh,
                       MPI_Op op );

  void singleFaceVectorFieldSymmetryBC( const int face,
                                        arrayView3d< real64 > const & vectorMultiField,
                                        arrayView3d< real64 > const & dVectorMultiField,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                        Group & nodeSets );

  void enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         Group & nodeSets );

  void applyEssentialBCs( const real64 dt,
                          const real64 time_n,
                          NodeManager & nodeManager );

  void applySuperimposedVelocityGradient( const real64 dt,
                                          ParticleManager & particleManager,
                                          SpatialPartition & partition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 bruteForceNodalAreaIntegration( real64 const (&hEl)[3],
                                         integer const & numSurfaceIntegrationPoints,
                                         real64 const & L,
                                         real64 const & dA,
                                         arraySlice1d< real64 const > const surfaceNormal,
                                         arraySlice1d< real64 const > const surfacePosition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 meshNodalAreaIntegration( real64 const (&hEl)[3],
                                  integer const & numSurfaceIntegrationPoints,
                                  localIndex const & numNeighbors,
                                  arraySlice1d< localIndex const > const neighborRegions,
                                  arraySlice1d< localIndex const > const neighborSubRegions,
                                  arraySlice1d< localIndex const > const neighborIndices,
                                  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleCenter,
                                  ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleRVectors,
                                  arraySlice1d< real64 const > const surfaceNormal,
                                  arraySlice1d< real64 const > const surfacePosition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 convexHullAreaIntegration( real64 const (& hEl)[3],
                                                     real64 const (& nodePosition)[3],
                                                     integer const & maxNodalNeighbors,
                                                     integer const & numSurfaceIntegrationPoints,
                                                     localIndex const & numNeighbors,
                                                     arraySlice1d< localIndex const > const neighborRegions,
                                                     arraySlice1d< localIndex const > const neighborSubRegions,
                                                     arraySlice1d< localIndex const > const neighborIndices,
                                                     ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleCenter,
                                                     ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleRVectors,
                                                     arraySlice1d< real64 const > const surfaceNormal,
                                                     arraySlice1d< real64 const > const surfacePosition );

  void computeNodalAreas( NodeManager & nodeManager );

  void normalizeGridSurfaceNormalsAndPositions( NodeManager & nodeManager );

  // void normalizeGridSurfacePositions( NodeManager & nodeManager );

  void computeGridSurfaceNormalWeights( ParticleManager & particleManager,
                                        NodeManager & nodeManager );

  void computeContactForces( real64 const dt,
                             ParticleManager & particleManager,
                             NodeManager & nodeManager );

  void initializeFrictionCoefficients();

  void lookupFrictionCoefficient( int const a,
                                  int const b,
                                  real64 & frictionCoefficient );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computePairwiseNodalContactForce( ContactGapCorrectionOption const & contactGapCorrection,
                                         OverlapCorrectionOption const & overlapCorrection,
                                         real64 const (&hEl) [3],
                                         int const & planeStrain,
                                         real64 const & smallMass,
                                         int const & useSurfacePositionForContact,
                                         int const & useCohesiveTangentialForces,
                                         int & separable,
                                         real64 const & dt,
                                         real64 const & frictionCoefficient,
                                         real64 ( &nAB )[3],
                                         real64 const & mA,
                                         real64 const & mB,
                                         real64 const & VA,
                                         real64 const & VB,
                                         arraySlice1d< real64 const > const vA,
                                         arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                         arraySlice1d< real64 const > const qA,
                                         arraySlice1d< real64 const > const qB,
                                         real64 const spmA,
                                         real64 const spmB,
                                         arraySlice1d< real64 const > const sA,
                                         arraySlice1d< real64 const > const sB,
                                         arraySlice1d< real64 const > const xA, // Position of field A
                                         arraySlice1d< real64 const > const xB, // Position of field B
                                         arraySlice1d< real64 const > const centerOfVolumeA,
                                         arraySlice1d< real64 const > const centerOfVolumeB,
                                         arraySlice1d< real64 > const fA,
                                         arraySlice1d< real64 > const fB );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeOrthonormalBasis( const real64 * e1,  // input "normal" unit vector.
                                real64 * e2,        // output "tangential" unit vector.
                                real64 * e3 );      // output "tangential" unit vector.

  void setGridFieldLabels( NodeManager & nodeManager );

  void solverProfiling( std::string label );

  real64 computeNeighborList( ParticleManager & particleManager );

  void generateNodalNeighborList( ParticleManager & particleManager,
                                  NodeManager & nodeManager );

  void optimizeBinSort( ParticleManager & particleManager );

  void particleColorSort( ParticleManager & particleManager );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernel( real64 const & r ); // distance from particle to query point

  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 inverseKernel( real64 const & d );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernelDevice( real64 const & r,
                       real64 const & neighborRadius,
                       int const & planeStrain ); // distance from particle to query point

  void kernelGradient( arraySlice1d< real64 const > const & x,  // query point
                       arraySlice1d< real64 const > const & xp, // particle location
                       real64 const & r,                        // distance from particle to query point
                       real64 ( &result )[3] );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void kernelGradientDevice( arraySlice1d< real64 const > const & x,  // query point
                             real64 const (&xp)[3],  //arraySlice1d< real64 const > const & xp, // particle location
                             real64 const & r,                        // distance from particle to query point
                             real64 const & neighborRadius,
                             int const & planeStrain,
                             real64 ( &result )[3] );

  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 computeKernelField( arraySlice1d< real64 const > const x,    // query point
                             localIndex const numNeighbors,
                             arrayView2d< real64 const > const & xp,    // List of neighbor particle locations.
                             arrayView1d< real64 const > const & Vp,    // List of neighbor particle volumes.
                             arrayView1d< real64 const > const & fp );  // scalar field values (e.g. damage) at neighbor particles

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeKernelFieldDevice( arraySlice1d< real64 const > const x,  // query point
                                   real64 const & neighborRadius,
                                   int const & planeStrain,
                                   localIndex const & numNeighbors,
                                   arraySlice1d< localIndex const > const regionIndices,
                                   arraySlice1d< localIndex const > const subRegionIndices,
                                   arraySlice1d< localIndex const > const particleIndices,
                                   ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                   ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                   const std::function< real64( localIndex, localIndex, localIndex ) > damage );

  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelFieldGradient( arraySlice1d< real64 const > const x, // query point
                                   localIndex const numNeighbors,
                                   arrayView2d< real64 const > const & xp,           // List of neighbor particle locations.
                                   arrayView1d< real64 const > const & Vp,           // List of neighbor particle volumes.
                                   arrayView1d< real64 const > const & fp,           // scalar field values (e.g. damage) at neighbor
                                                                                     // particles
                                   arraySlice1d< real64 > const result );

  void computeKernelFieldGradient( arraySlice1d< real64 const > const x,       // query point
                                   std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                   std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                   std::vector< real64 > & fp,                 // scalar field values (e.g. damage) at neighbor particles
                                   real64 ( &result )[3] );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldGradientDevice( arraySlice1d< real64 const > const x, // query point
                                         real64 const & neighborRadius,
                                         int const & planeStrain,
                                         localIndex const & numNeighbors,
                                         arraySlice1d< localIndex const > const regionIndices,
                                         arraySlice1d< localIndex const > const subRegionIndices,
                                         arraySlice1d< localIndex const > const particleIndices,
                                         ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                         ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                         const std::function< real64( localIndex, localIndex, localIndex ) > damage,
                                         real64 ( &result )[3] );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldDivergenceDevice( arraySlice1d< real64 const > const x, // query point
                                           real64 const & neighborRadius,
                                           int const & planeStrain,
                                           localIndex const & numNeighbors,
                                           arraySlice1d< localIndex const > const regionIndices,
                                           arraySlice1d< localIndex const > const subRegionIndices,
                                           arraySlice1d< localIndex const > const particleIndices,
                                           ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particleVectorView,
                                           real64 & result );

  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelVectorGradient( arraySlice1d< real64 const > const x,       // query point
                                    localIndex const numNeighbors,
                                    arrayView2d< real64 const > const & xp,  // List of neighbor particle locations.
                                    arrayView1d< real64 const > const & Vp,                 // List of neighbor particle volumes.
                                    arrayView2d< real64 const > const & fp,  // vector field values (e.g. velocity) at neighbor particles
                                    arraySlice2d< real64 > const result );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelVectorGradientDevice( arraySlice1d< real64 const > const x,       // query point
                                          real64 const & neighborRadius,
                                          int const & planeStrain,
                                          localIndex const numNeighbors,
                                          arraySlice1d< localIndex const > const regionIndices,
                                          arraySlice1d< localIndex const > const subRegionIndices,
                                          arraySlice1d< localIndex const > const particleIndices,
                                          ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                          ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                          ParticleManager::ParticleView< arrayView2d< real64 const > > particleDisplacementView,
                                          arraySlice2d< real64 > const result );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  bool markSurfaceAsDamage( int const & surfaceFlag );

  void computeDamageFieldGradient( ParticleManager & particleManager );

  void computeDamageHessian( ParticleManager & particleManager );

  void tagBinderCZSurfaces( ParticleManager & particleManager,
                            NodeManager & nodeManager );

  void czSurfaceFlagUpdate( ParticleManager & particleManager );

  void computeDistanceToCrackTip( ParticleManager & particleManager );

  void estimateParticleSurfacePosition( ParticleManager & particleManager );

  void updateSurfaceFlagOverload( ParticleManager & particleManager );

  void projectDamageFieldGradientToGrid( DomainPartition & domain,
                                         ParticleManager & particleManager,
                                         NodeManager & nodeManager,
                                         MeshLevel & mesh );

  void computeParticleFieldMappings( ParticleManager & particleManager,
                                     NodeManager & nodeManager );

  void updateDeformationGradient( real64 dt,
                                  ParticleManager & particleManager );

  void stressControl( const real64 dt,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  void initializeParticleFields( ParticleManager & particleManager );

  void initializeConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateStress( real64 dt,
                     ParticleManager & particleManager );

  void particleKinematicUpdate( const real64 dt,
                                ParticleManager & particleManager );

  void computeAndWriteBoxAverage( const real64 dt,
                                  const real64 time_n,
                                  ParticleManager & particleManager );

  void writeParticleData( const real64 time_n,
                          ParticleManager & particleManager );

  void computeBoxMetrics( ParticleManager & particleManager,
                          arrayView1d< real64 > boxStress,
                          real64 & boxMaterialVolume );

  void initializeGridFields( NodeManager & nodeManager );

  void boundaryConditionUpdate( real64 dt, real64 time_n );

  void projectParticleSurfaceNormalsToGrid( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            MeshLevel & mesh );

  void initializeCohesiveReferenceConfiguration( DomainPartition & domain,
                                                 ParticleManager & particleManager,
                                                 NodeManager & nodeManager,
                                                 MeshLevel & mesh );

  bool interiorToParticleProjectedArea( ParticleManager & particleManager,
                                        globalIndex const GEOS_UNUSED_PARAM( gridIndex ),
                                        int const gridFieldIndex,
                                        real64 const (&gridSurfaceNormal)[3],
                                        real64 const (&gridSurfacePoint)[3] );

  void projectToPlane( real64 const (&vector)[3],
                       real64 const (&normal)[3],
                       real64 ( &projection )[3] );

  void enforceCohesiveLaw( real64 dt,
                           ParticleManager & particleManager,
                           NodeManager & nodeManager );

  GEOS_HOST_DEVICE
  real64 computeDistanceToParticleSurface( real64 ( &normal )[3],
                                           arraySlice2d< real64 const > const rVectors );

  void updateGridBackgroundStress( NodeManager & nodeManager );

  void particleToGrid( real64 const time_n,
                       integer const cycleNumber,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  void particleToGrid_reduction( real64 const time_n,
                                 integer const cycleNumber,
                                 ParticleManager & particleManager,
                                 NodeManager & nodeManager );

  void particleToGrid_noAtomics( real64 const time_n,
                                 integer const cycleNumber,
                                 ParticleManager & particleManager,
                                 NodeManager & nodeManager );

  void particleToGrid_randomMix( real64 const time_n,
                                 integer const cycleNumber,
                                 ParticleManager & particleManager,
                                 NodeManager & nodeManager );

  void particleToGrid_minimalAtomics( real64 const time_n,
                                      integer const cycleNumber,
                                      ParticleManager & particleManager,
                                      NodeManager & nodeManager );

  void particleToGrid_colors( real64 const time_n,
                              integer const cycleNumber,
                              ParticleManager & particleManager,
                              NodeManager & nodeManager );

  void computeSPHSurfaceCurvature( ParticleManager & particleManager );

  void gridTrialUpdate( real64 dt,
                        NodeManager & nodeManager );

  void enforceContact( real64 dt,
                       //  DomainPartition & domain,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager
                       //  MeshLevel & mesh
                       );

  void computeParticleSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                  NodeManager & nodeManager );

  void computePairwiseLogisticRegressionSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                                    NodeManager & nodeManager );

  void computeGridSurfacePositionFromVolume( NodeManager & nodeManager );

  real64 computeVolumeFromSurfacePosition( real64 const (&n)[3],
                                           real64 const offset,
                                           real64 const (&hEl) [3] );

  real64 computeMaximumSurfacePositionOffset( real64 ( &n )[3],
                                              real64 const (&hEl)[3] );

  void logisticRegression( int const & planeStrain,
                           integer const & numContactGroups,
                           integer const & damageFieldPartitioning,
                           integer const & maxLRIterations,
                           real64 const & LRtolerance,
                           real64 const (&hEl)[3],
                           localIndex const & fieldA,
                           localIndex const & fieldB,
                           localIndex const & numNeighboringParticles,
                           arraySlice1d< localIndex const > const neighborRegions,
                           arraySlice1d< localIndex const > const neighborSubRegions,
                           arraySlice1d< localIndex const > const neighborIndices,
                           ParticleManager::ParticleViewConst< arrayView1d< localIndex const > > particleGroup,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleDamageGradient,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleSurfaceNormal,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePosition,
                           arraySlice1d< real64 const > const gridPosition,
                           arraySlice1d< real64 const > const gridDamageGradient,
                           real64 ( &normal0 )[3],
                           real64 ( &normal )[3],
                           real64 ( &surfacePosition )[3] );

  void mapSurfaceNormalsAndPositionsToParticles( ParticleManager & particleManager,
                                                 NodeManager & nodeManager );

  void interpolateTable( real64 const & x,
                         real64 const & dx,
                         arrayView2d< real64 const > const &  table,
                         arrayView1d< real64 > output,
                         arrayView1d< real64 > outputRate,
                         mpm::InterpolationOption const & interpolationType );

  void interpolateValueInRange( real64 const & x,
                                real64 const & xmin,
                                real64 const & xmax,
                                real64 const & ymin,
                                real64 const & ymax,
                                real64 & output,
                                mpm::InterpolationOption interpolationType );

  void interpolateFTable( real64 dt, real64 time_n );

  void interpolateStressTable( real64 dt, real64 time_n );

  void interpolateTemperatureTable( real64 dt, real64 time_n );

  void gridToParticle( real64 dt,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       DomainPartition & domain,
                       MeshLevel & mesh );

  void performFLIPUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager );

  void performPICUpdate( real64 dt,
                         ParticleManager & particleManager,
                         NodeManager & nodeManager );

  void performXPICUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  void performFMPMUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  void updateSolverDependencies( ParticleManager & particleManager );

  real64 getStableTimeStep( ParticleManager & particleManager );

  void deleteBadParticles( ParticleManager & particleManager );

  void printProfilingResults();

  void computeSurfaceFlags( ParticleManager & particleManager );

  void computeSurfaceNormals( ParticleManager & particleManager,
                              NodeManager & nodeManager );

  void computeSurfacePositions( ParticleManager & particleManager,
                                NodeManager & nodeManager );

  void computeSphF( ParticleManager & particleManager );

  // void directionalOverlapCorrection( real64 dt, ParticleManager & particleManager );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  int evaluateSeparabilityCriterion( int const & planeStrain,
                                     int const & numContactGroups,
                                     int const & treatFullyDamagedAsSingleField,
                                     real64 const & separabilityMinDamage,
                                     real64 const & thinFeatureDFGThreshold,
                                     real64 const & neighborRadius,
                                     localIndex const & A,
                                     localIndex const & B,
                                     real64 const & damageA,
                                     real64 const & damageB,
                                     real64 const & maxDamageA,
                                     real64 const & maxDamageB,
                                     arraySlice1d< real64 const > const damageGradient,
                                     arraySlice1d< real64 const > const xA,
                                     arraySlice1d< real64 const > const xB );

  void flagOutOfRangeParticles( ParticleManager & particleManager );

  void computeRVectors( ParticleManager & particleManager );

  void cpdiDomainScaling( ParticleManager & particleManager );

  void subdivideParticles( ParticleManager & particleManager );

  void resizeMappingArrays( ParticleManager & particleManager );

  void populateMappingArrays( ParticleManager & particleManager,
                              NodeManager & nodeManager ); //,
  //  SpatialPartition & partition  );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeSinglePointShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         arraySlice1d< real64 const > const particlePosition,
                                         arrayView3d< int const > const ijkMap,
                                         real64 const (&xLocalMin)[3],
                                         real64 const (&hEl)[3],
                                         arraySlice1d< int > const mappedNodes,
                                         arraySlice1d< real64 > const shapeFunctionValues,
                                         arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeCPDIShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                  arraySlice1d< real64 const > const particlePosition,
                                  arraySlice2d< real64 const > const particleRVectors,
                                  arrayView3d< int const > const ijkMap,
                                  real64 const (&xLocalMin)[3],
                                  real64 const (&hEl)[3],
                                  arraySlice1d< int > const mappedNodes,
                                  arraySlice1d< real64 > const shapeFunctionValues,
                                  arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void mapNodesAndComputeShapeFunctions( arrayView3d< localIndex const > const ijkMap,
                                         real64 const (&xLocalMin)[3],
                                         real64 const (&hEl)[3],
                                         ParticleType particleType,
                                         arraySlice1d< real64 const > const particlePosition,
                                         arraySlice2d< real64 const > const particleRVectors,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         localIndex * const mappedNode,
                                         real64 * const shapeFunctionValues,
                                         real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeSinglePointParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                 arraySlice1d< real64 const > const particlePosition,
                                                 arrayView3d< int const > const ijkMap,
                                                 real64 const (&xLocalMin)[3],
                                                 real64 const (&hEl)[3],
                                                 localIndex * const mappedNodes,
                                                 real64 * const shapeFunctionValues,
                                                 real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeCPDIParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                          arraySlice1d< real64 const > const particlePosition,
                                          arraySlice2d< real64 const > const particleRVectors,
                                          arrayView3d< int const > const ijkMap,
                                          real64 const (&xLocalMin)[3],
                                          real64 const (&hEl)[3],
                                          localIndex * const mappedNodes,
                                          real64 * const shapeFunctionValues,
                                          real64 shapeFunctionGradientValues[][3] );

  void computeBodyForce( ParticleManager & particleManager );

  void computeArtificialViscosity( ParticleManager & particleManager );

  void computeSPHJacobian( ParticleManager & particleManager );

  void overlapCorrection( real64 const dt,
                          ParticleManager & particleManager );

  void sphOverlapCorrection( real64 const dt,
                             ParticleManager & particleManager );

  void computeInternalEnergyAndTemperature( const real64 dt,
                                            ParticleManager & particleManager );

  void computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                             ParticleManager & particleManager );

  void correctGhostParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                            SpatialPartition & partition );

  void correctParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                       SpatialPartition & partition );

  void resetDeformationGradient( ParticleManager & particleManager );

  void unscaleCPDIVectors( ParticleManager & particleManager );

  void computeKineticEnergy( ParticleManager & particleManager );

  void computeXProfile( int const cycleNumber,
                        real64 const time,
                        real64 const dt,
                        NodeManager & nodeManager,
                        SpatialPartition & partition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void cofactor( real64 const (&F)[3][3],
                 real64 ( &Fc )[3][3] );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 Mod( real64 num, real64 denom );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int combinations( int n,
                    int k );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int factorial( int n );

protected:
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode );

  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode,
                                  xmlWrapper::xmlNodePos const & targetNodePos );

  virtual void postInputInitialization() override final;

  virtual void postRestartInitialization() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  stdVector< array2d< localIndex > > m_mappedNodes; // mappedNodes[subregion index][particle index][node index]. dims = {# of subregions,
                                                    // # of particles, # of nodes a particle on the subregion maps to}
  stdVector< array2d< integer > > m_mappedFields;
  stdVector< array2d< real64 > > m_shapeFunctionValues; // mappedNodes[subregion][particle][nodal shape function value]. dims = {# of
                                                        // subregions, # of particles, # of nodes a particle on the subregion maps to}
  stdVector< array3d< real64 > > m_shapeFunctionGradientValues; // mappedNodes[subregion][particle][nodal shape function gradient
                                                                // value][direction]. dims = {# of subregions, # of particles, # of nodes
                                                                // a particle on the subregion maps to, 3}

  // Domain options
  int m_numDims;
  int m_planeStrain;
  array1d< int > m_nEl;                   // Number of elements in each grid direction including buffer and ghost cells
  array1d< real64 > m_hEl;                // Grid spacing in x-y-z
  array1d< real64 > m_xLocalMin;          // Minimum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMax;          // Maximum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMinNoGhost;   // Minimum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xLocalMaxNoGhost;   // Maximum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xGlobalMin;         // Minimum global grid coordinate excluding buffer nodes
  array1d< real64 > m_xGlobalMax;         // Maximum global grid coordinate excluding buffer nodes
  array1d< real64 > m_partitionExtent;    // Length of each edge of partition including buffer and ghost cells
  array1d< real64 > m_domainExtent;       // Length of each edge of global domain excluding buffer cells
  array3d< localIndex > m_ijkMap;        // Map from indices in each spatial dimension to local node ID
  localIndex m_numberOfSubRegions;
  real64 m_smallMass;

  // Debugging / Profiling options
  int m_solverProfiling;
  int m_logStartCycle;
  stdVector< real64 > m_profilingTimes;
  stdVector< std::string > m_profilingLabels;

  // Plotting options
  int m_plotGridFields;
  string_array m_plottableFields;
  SortedArray< string > m_plottableFieldsSorted;
  int m_plotUnscaledParticles;

  // Output options
  int m_boxAverageHistory;
  int m_boxAverageResizeWithDomain;  // 0 use constant box_average domain, 1: resize with domainL
  real64 m_boxAverageWriteInterval;
  real64 m_nextBoxAverageWriteTime;
  array1d< real64 > m_boxAverageMin;
  array1d< real64 > m_boxAverageMax;

  int m_reactionHistory;
  real64 m_reactionWriteInterval;
  real64 m_nextReactionWriteTime;

  int m_writeParticleData;
  real64 m_particleDataWriteInterval;
  real64 m_nextParticleDataWriteTime;

  int m_computeXProfile;
  real64 m_xProfileWriteInterval;
  real64 m_nextXProfileWriteTime;
  real64 m_xProfileVx0;

  // Integration options
  TimeIntegrationOption m_timeIntegrationOption;
  UpdateMethodOption m_updateMethod;
  int m_updateOrder;
  int m_FSubcycles;
  int m_LBar;
  real64 m_LBarScale;
  int m_exactJIntegration;
  int m_useAPIC;
  int m_useInteralForceAsFaceReaction;

  // Boundary condition options
  int m_prescribedBcTable;
  array2d< real64 > m_bcTable;
  array1d< int > m_boundaryConditionTypes; // TODO: Surely there's a way to have just one variable here
  array1d< real64 > m_boundaryFaceCoefficientsOfRestitution;
  array1d< real64 > m_boundaryFaceFrictionCoefficients; // Ignored unless face has boundary condition type 3

  int m_prescribedFTable;
  int m_prescribedBoundaryFTable;
  mpm::InterpolationOption m_fTableInterpType;
  array2d< real64 > m_fTable;
  array1d< real64 > m_domainF;
  array1d< real64 > m_domainL;

  array1d< int > m_enablePrescribedBoundaryTransverseVelocities;
  array2d< real64 > m_prescribedBoundaryTransverseVelocities; // 2 in-plane directions * 6 faces

  array1d< real64 > m_globalFaceReactions;

  array1d< int > m_stressControl;
  mpm::InterpolationOption m_stressTableInterpType;
  array2d< real64 > m_stressTable;
  real64 m_stressControlKp;
  real64 m_stressControlKi;
  real64 m_stressControlKd;
  array1d< real64 > m_stressControlLastError;
  array1d< real64 > m_stressControlITerm;
  array1d< real64 > m_domainStress;

  array1d< real64 > m_bodyForce;

  int m_enableSurfaceTension;
  real64 m_surfaceTensionCoefficient;

  // Borehole options - borehole fluid pressure and radius used in the boreholePressure event.
  int m_enableBoreholePressure;
  array1d< real64 > m_boreholeStress;
  real64 m_boreholeRadius;

  int m_enableConfiningPressure;
  array1d< real64 > m_confiningStress;
  array1d< real64 > m_confiningPressureBoxMin;
  array1d< real64 > m_confiningPressureBoxMax;

  // Temperature options
  real64 m_domainTemperature;
  real64 m_domainTemperatureRate;

  int m_shockHeating;
  int m_computeInternalEnergyAndTemperature;
  int m_useArtificialViscosity;
  real64 m_artificialViscosityQ0;
  real64 m_artificialViscosityQ1;

  real64 m_damageHessianSurfaceThreshold;
  int m_computeCZInterfacesFromDamage;
  int m_checkForBinder;

  // DFG options
  int m_damageFieldPartitioning;
  int m_numContactGroups;
  int m_numContactFlags;
  int m_numVelocityFields;
  real64 m_separabilityMinDamage;
  int m_treatFullyDamagedAsSingleField;
  real64 m_thinFeatureDFGThreshold;
  int m_useDamageAsSurfaceFlag;

  array1d< int > m_groupNormalPriority;

  // Contact options
  int m_enableContact;
  int m_hasContact;
  ContactNormalTypeOption m_contactNormalType;
  ContactGapCorrectionOption m_contactGapCorrection;
  int m_useSurfacePositionForContact;
  real64 m_explicitSurfaceNormalInfluence;
  int m_computeParticleSurfaceNormalsAndPositions;
  real64 m_frictionCoefficient;
  array2d< real64 > m_frictionCoefficientTable;
  NormalsAndPositionsMethodOption m_normalAndPositionMethod;
  int m_overwriteExistingNormalsAndPositions;
  real64 m_contactNormalExponent;
  int m_maxLRIterations;
  real64 m_LRtolerance;
  int m_surfaceDetection;

  // Performance options
  GPUSchemeOption m_gpuScheme;

  // Cohesive zone options
  AreaIntegrationOption m_areaIntegrationMethod;
  real64 m_numSurfaceIntegrationPoints;
  int m_preventCZInterpenetration;
  real64 m_totalBinderVolume;
  real64 m_polymerCZThickness;

  // Neighborlist options
  int m_needsNeighborList;
  real64 m_neighborRadius;
  int m_binSizeMultiplier;
  int m_needsNodalNeighborList;
  OrderedVariableToManyParticleRelation m_nodalNeighborList;
  int m_maxNodalNeighbors;

  // Conditioning options
  int m_cpdiDomainScaling;
  int m_subdivideParticles; // Gas particles larger than a grid cell are subdivided
  real64 m_maxParticleVelocity;
  real64 m_maxParticleVelocitySquared;
  real64 m_minParticleJacobian;
  real64 m_maxParticleJacobian;
  int m_disableSurfaceNormalsAndPositionsOnCPDIScaling; // Turns off surface normals and positions for highly deformed particles
  int m_disableSurfaceNormalsAndPositionsOnDamage; // Turns off surface normals and positions for highly damaged particles
  real64 m_surfaceNormalAndPositionDamageThreshold;
  int m_resetDefGradForFullyDamagedParticles;
  int m_resetDefGradForScaledSurfaceParticles;
  int m_useReferenceVectorsForParticleUpdate;

  // Overlap correction options
  OverlapCorrectionOption m_overlapCorrection;
  real64 m_overlapThreshold1;
  real64 m_overlapThreshold2;
  int m_computeSPHJacobian;
  int m_directionalOverlapCorrection;

  // Crack Tip options - parameters for crack-tip detection used for stress concentration factor
  int m_useCrackTipDetection;
  real64 m_crackTipDetectionThreshold;

  // Event options
  int m_useEvents;                   // Events flag
  int m_surfaceHealing;

  // Method of manufactured solutions (MMS) options
  int m_generalizedVortexMMS;

  // Misc
  int m_computeNodalArea;
  int m_useNodePosForArea;

private:
  struct BinKey
  {
    localIndex regionIndex;
    localIndex subRegionIndex;
    localIndex binIndex;

    bool operator==( BinKey const & other ) const
    {
      return (regionIndex == other.regionIndex && subRegionIndex == other.subRegionIndex && binIndex == other.binIndex);
    }
  };

  struct BinKeyHash
  {
    std::size_t operator()( BinKey const & k ) const
    {
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      return ((std::hash< localIndex >()( k.regionIndex )
               ^ (std::hash< localIndex >()( k.subRegionIndex ) << 1)) >> 1)
             ^ (std::hash< localIndex >()( k.binIndex ) << 1);
    }
  };

  void setParticlesConstitutiveNames( ParticleSubRegionBase & subRegion ) const;
};

ENUM_STRINGS( SolidMechanicsMPM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

ENUM_STRINGS( SolidMechanicsMPM::UpdateMethodOption,
              "FLIP",
              "PIC",
              "XPIC",
              "FMPM" );

ENUM_STRINGS( SolidMechanicsMPM::BoundaryConditionOption,
              "Outflow",
              "Symmetry",
              "Moving",
              "Contact" );

ENUM_STRINGS( SolidMechanicsMPM::ContactNormalTypeOption,
              "Difference",
              "MassWeighted",
              "LargerMass",
              "Mixed",
              "Aligned",
              "LogisticRegression" );

ENUM_STRINGS( SolidMechanicsMPM::ContactGapCorrectionOption,
              "Simple",
              "Implicit",
              "Softened" );

ENUM_STRINGS( SolidMechanicsMPM::AreaIntegrationOption,
              "BruteForce",
              "Mesh" );

ENUM_STRINGS( SolidMechanicsMPM::OverlapCorrectionOption,
              "Off",
              "NormalForce",
              "SPH" );

ENUM_STRINGS( SolidMechanicsMPM::CohesiveLawOption,
              "Uncoupled",
              "NeedlemanXu",
              "Polymer" );

ENUM_STRINGS( SolidMechanicsMPM::GPUSchemeOption,
              "Atomics",
              "NoAtomics",
              "RandomMix",
              "MinimalAtomics",
              "Reduction",
              "Colors" );

ENUM_STRINGS( SolidMechanicsMPM:: NormalsAndPositionsMethodOption,
              "LogisticRegression",
              "DFGAndVolumeIntegration" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMPM_HPP_ */
