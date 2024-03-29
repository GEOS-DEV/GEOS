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
 * @file SolidMechanicsMPM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "common/TimingMacros.hpp"
#include "kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "kernels/ExplicitMPM.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "MPMSolverFields.hpp"

#include "events/mpmEvents/MPMEventManager.hpp"

namespace geos
{

class SpatialPartition;


/**
 * @class SolidMechanicsMPM
 *
 * This class implements a material point method solution to the equations of motion.
 */
class SolidMechanicsMPM : public SolverBase
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
    OUTFLOW,    //!<Outflow
    SYMMETRY    //!<Symmetry
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
    Cohesive    
  };

  /**
   * @enum InterpolationOption
   * 
   * The options for interpolating tables
  */
  enum struct InterpolationOption : integer
  {
    Linear,
    Cosine,
    Smoothstep
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
    Mixed
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
    PPR,
    NeedlemanXu
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
  virtual ~SolidMechanicsMPM() override {};

  /**
   * @return The string that may be used to generate a new instance from the SolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanics_MPM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  void updateIntrinsicNodalData( DomainPartition * const domain );


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
                                 arrayView1d< string const > const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  /**
   * Applies displacement boundary conditions to the system for implicit time integration
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * forceExternalString() { return "externalForce"; }
    static constexpr char const * forceInternalString() { return "internalForce"; }
    static constexpr char const * displacementString() { return "displacement"; }
    static constexpr char const * particleSurfaceNormalString() { return "particleSurfaceNormal"; }
    static constexpr char const * cohesiveTractionString() { return "cohesiveTraction"; }

    static constexpr char const * massString() { return "mass"; }
    static constexpr char const * materialVolumeString() { return "materialVolume"; }
    static constexpr char const * velocityString() { return "velocity"; }
    static constexpr char const * dVelocityString() { return "dVelocity"; }
    static constexpr char const * momentumString() { return "momentum"; }
    static constexpr char const * accelerationString() { return "acceleration"; }
    static constexpr char const * forceContactString() { return "contactForce"; }
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * damageGradientString() { return "damageGradient"; }
    static constexpr char const * maxDamageString() { return "maxDamage"; }
    static constexpr char const * surfaceNormalString() { return "surfaceNormal"; }

    static constexpr char const * surfacePositionString() { return "surfacePosition"; }

    static constexpr char const * materialPositionString() { return "materialPosition"; }
    static constexpr char const * normalStressString() { return "normalStress"; }
    static constexpr char const * massWeightedDamageString() { return "massWeightedDamage"; }
    static constexpr char const * cohesiveNodeString() { return "cohesiveNode"; }
    static constexpr char const * referenceAreaVectorString() { return "referenceAreaVector"; }
    static constexpr char const * referenceSurfacePositionString() { return "referenceSurfacePosition"; }
    static constexpr char const * referenceMaterialVolumeString() { return "referenceMaterialVolume"; }

    static constexpr char const * surfaceMassString() { return "surfaceMass"; }
    static constexpr char const * surfaceFieldMassString() { return "surfaceFieldMass"; }
    static constexpr char const * explicitSurfaceNormalString() { return "explicitSurfaceNormal"; }
    static constexpr char const * maxMappedParticleIDString() { return "maxMappedParticleIDS"; }
    static constexpr char const * principalExplicitSurfaceNormalString() { return "principalExplicitSurfaceNormal"; }
    static constexpr char const * cohesiveFieldFlagString() { return "cohesiveFieldFlag"; }
    static constexpr char const * cohesiveForceString() { return "cohesiveForce"; }
    
    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }

    static constexpr char const * vPlusString() { return "vPlus"; }
    static constexpr char const * dVPlusString() { return "dVPlus"; }

    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  /// Child group viewKeys
  struct groupKeysStruct
  {
    dataRepository::GroupKey mpmEventManager = { "MPMEvents" }; ///< MPM Events key
  } groupKeys; ///< Child group viewKeys

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

  void performMaterialSwap( ParticleManager & particleManager,
                            string sourceRegionName,
                            string destinationRegionName );

  void resizeGrid( SpatialPartition & partition,
                   NodeManager & nodeManager,
                   real64 const dt );

  void syncGridFields( std::vector< std::string > const & fieldNames,
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

  void computeGridSurfaceNormals( ParticleManager & particleManager,
                                  NodeManager & nodeManager );

  void computeGridSurfacePositions( ParticleManager & particleManager,
                                    NodeManager & nodeManager );

  void normalizeGridSurfaceNormals( NodeManager & nodeManager );

  void normalizeGridSurfacePositions( NodeManager & nodeManager );

  void computeContactForces( real64 const dt,
                             NodeManager & nodeManager,
                             arrayView2d< real64 const > const & gridMass,
                             arrayView2d< real64 const > const & gridMaterialVolume,
                             arrayView2d< real64 const > const & gridDamage,
                             arrayView2d< real64 const > const & gridMaxDamage,
                             arrayView3d< real64 const > const & gridVelocity,
                             arrayView3d< real64 const > const & gridMomentum,
                             arrayView3d< real64 const > const & gridSurfaceNormal,
                             arrayView3d< real64 const > const & gridSurfacePosition,
                             arrayView3d< real64 const > const & gridMaterialPosition,
                             arrayView2d< int const > const & gridCohesiveFieldFlag,
                             arrayView3d< real64 > const & gridContactForce );

  void initializeFrictionCoefficients();

  void lookupFrictionCoefficient( int const a,
                                  int const b,
                                  real64 & frictionCoefficient );

  void computePairwiseNodalContactForce( int const & separable,
                                         int useCohesiveTangentialForces,
                                         real64 const & dt,
                                         real64 const & frictionCoefficient,
                                         real64 const & mA,
                                         real64 const & mB,
                                         real64 const & VA,
                                         real64 const & VB,
                                         arraySlice1d< real64 const > const vA,
                                         arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                         arraySlice1d< real64 const > const qA,
                                         arraySlice1d< real64 const > const qB,
                                         arraySlice1d< real64 const > const nA,
                                         arraySlice1d< real64 const > const nB,
                                         arraySlice1d< real64 const > const sA,
                                         arraySlice1d< real64 const > const sB, 
                                         arraySlice1d< real64 const > const xA, // Position of field A
                                         arraySlice1d< real64 const > const xB, // Position of field B
                                         arraySlice1d< real64 > const fA,
                                         arraySlice1d< real64 > const fB );

  void computeOrthonormalBasis( const real64 * e1,  // input "normal" unit vector.
                                real64 * e2,        // output "tangential" unit vector.
                                real64 * e3 );      // output "tangential" unit vector.

  void setGridFieldLabels( NodeManager & nodeManager );

  void solverProfiling( std::string label );

  void solverProfilingIf( std::string label, bool condition );

  real64 computeNeighborList( ParticleManager & particleManager );

  void optimizeBinSort( ParticleManager & particleManager );

  real64 kernel( real64 const & r ); // distance from particle to query point

  void kernelGradient( arraySlice1d< real64 const > const x,  // query point
                       std::vector< real64 > & xp,            // particle location
                       real64 const & r,                      // distance from particle to query point
                       real64 * result );

  real64 computeKernelField( arraySlice1d< real64 const > const x,    // query point
                             arrayView2d< real64 const > const xp,    // List of neighbor particle locations.
                             arrayView1d< real64 const > const Vp,    // List of neighbor particle volumes.
                             arrayView1d< real64 const > const fp );  // scalar field values (e.g. damage) at neighbor particles

  void computeKernelFieldGradient( arraySlice1d< real64 const > const x,       // query point
                                   std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                   std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                   std::vector< real64 > & fp,                 // scalar field values (e.g. damage) at neighbor particles
                                   arraySlice1d< real64 > const result );

  void computeKernelVectorGradient( arraySlice1d< real64 const > const x,       // query point
                                    std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                    std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                    std::vector< std::vector< real64 > > & fp,  // vector field values (e.g. velocity) at neighbor particles
                                    arraySlice2d< real64 > const result );

  void computeDamageFieldGradient( ParticleManager & particleManager );

  void updateSurfaceFlagOverload( ParticleManager & particleManager );

  void projectDamageFieldGradientToGrid( DomainPartition & domain,
                                         ParticleManager & particleManager,
                                         NodeManager & nodeManager,
                                         MeshLevel & mesh );

  void updateDeformationGradient( real64 dt,
                                  ParticleManager & particleManager );

  void stressControl( const real64 dt,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  void initializeConstitutiveModelDependencies( ParticleManager & particleManager);

  void updateConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateStress( real64 dt,
                     ParticleManager & particleManager );

  void particleKinematicUpdate( ParticleManager & particleManager );

  void computeAndWriteBoxAverage( const real64 dt,
                                  const real64 time_n,
                                  ParticleManager & particleManager );

  void computeBoxMetrics( ParticleManager & particleManager,
                          arrayView1d< real64 > boxStress,
                          real64 & boxMaterialVolume );

  void initializeGridFields( NodeManager & nodeManager );

  void boundaryConditionUpdate( real64 dt, real64 time_n );

  void projectParticleSurfaceNormalsToGrid( DomainPartition & domain,
                                            ParticleManager& particleManager,
                                            NodeManager & nodeManager,
                                            MeshLevel & mesh  );

  void initializeCohesiveReferenceConfiguration( DomainPartition & domain,
                                                 ParticleManager& particleManager,
                                                 NodeManager & nodeManager,
                                                 MeshLevel & mesh );

  bool interiorToParticleProjectedArea( ParticleManager & particleManager,
                                        globalIndex const GEOS_UNUSED_PARAM( gridIndex ),
                                        int const gridFieldIndex,
                                        real64 const (& gridSurfaceNormal)[3],
                                        real64 const (& gridSurfacePoint)[3] );

  void projectToPlane( real64 const (& vector)[3],
                       real64 const (& normal)[3],
                       real64 (& projection)[3] );

  void enforceCohesiveLaw(  ParticleManager & particleManager,
                            NodeManager & nodeManager );

  void computeDistanceToParticleSurface( real64 (& normal)[3],
                                         arraySlice2d< real64 const > const rVectors,
                                         real64 distanceToSurface );

  void computeCohesiveTraction( int g,
                                int a,
                                int b, 
                                real64 mA,
                                real64 mB,
                                arraySlice1d< real64 const > const dA,
                                arraySlice1d< real64 const > const dB,
                                real64 const (& aA )[3],
                                real64 const (& aB )[3],
                                arraySlice1d< real64 const > const nA,
                                arraySlice1d< real64 const > const nB, 
                                arraySlice1d< real64 > const tA,
                                arraySlice1d< real64 > const tB );

void uncoupledCohesiveLaw( real64 normalDisplacement,
                           real64 tangentialDisplacement,
                           real64 & normalStress,
                           real64 & shearStress );

void needlemanXuCohesiveLaw( real64 normalDisplacement,
                             real64 tangentialDisplacement,
                             real64 & normalStress,
                             real64 & shearStress );

  void particleToGrid( real64 const time_n,
                       integer const cycleNumber,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  void gridTrialUpdate( real64 dt,
                        NodeManager & nodeManager );

  void enforceContact( real64 dt,
                       DomainPartition & domain,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       MeshLevel & mesh );

void interpolateTable( real64 x, 
                       real64 dx,
                       array2d< real64 > table,
                       arrayView1d< real64 > output,
                       SolidMechanicsMPM::InterpolationOption interpolationType );

  void interpolateFTable( real64 dt, real64 time_n );

  void interpolateStressTable( real64 dt, real64 time_n );

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

  void computeSphF( ParticleManager & particleManager );

  // void directionalOverlapCorrection( real64 dt, ParticleManager & particleManager );

  int evaluateSeparabilityCriterion( localIndex const & A,
                                     localIndex const & B,
                                     real64 const & damageA,
                                     real64 const & damageB,
                                     real64 const & maxDamageA,
                                     real64 const & maxDamageB,
                                     arraySlice1d< real64 const > const xA,
                                     arraySlice1d< real64 const > const xB );

  void flagOutOfRangeParticles( ParticleManager & particleManager );

  void computeRVectors( ParticleManager & particleManager );

  void cpdiDomainScaling( ParticleManager & particleManager );

  void resizeMappingArrays( ParticleManager & particleManager );

  void populateMappingArrays( ParticleManager & particleManager,
                              NodeManager & nodeManager ); //,
                                              //  SpatialPartition & partition  );

GEOS_FORCE_INLINE
GEOS_HOST_DEVICE
void computeSinglePointShapeFunctions( arrayView2d< real64 const > const gridPosition,
                                       arraySlice1d< real64 const > const particlePosition,
                                       arrayView3d< int const > const ijkMap,
                                       real64 const (& xLocalMin)[3],
                                       real64 const (&hEl)[3],
                                       arraySlice1d< int > const mappedNodes,
                                       arraySlice1d< real64 > const shapeFunctionValues,
                                       arraySlice2d< real64 > const shapeFunctionGradientValues );

GEOS_FORCE_INLINE
GEOS_HOST_DEVICE
void computeCPDIShapeFunctions( arrayView2d< real64 const > const gridPosition,
                                arraySlice1d< real64 const > const particlePosition,
                                arraySlice2d< real64 const > const particleRVectors,
                                arrayView3d< int const > const ijkMap,
                                real64 const (& xLocalMin)[3],
                                real64 const (&hEl)[3],
                                arraySlice1d< int > const mappedNodes,
                                arraySlice1d< real64 > const shapeFunctionValues,
                                arraySlice2d< real64 > const shapeFunctionGradientValues );

  void computeBodyForce( ParticleManager & particleManager );
  
  void computeWavespeed( ParticleManager & particleManager );

  void computeArtificialViscosity( ParticleManager & particleManager );

  void computeSPHJacobian( ParticleManager & particleManager );

  void overlapCorrection( real64 const dt ,
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

  void cofactor( real64 const (& F)[3][3],
                 real64 (& Fc)[3][3] );

  real64 Mod(real64 num, real64 denom);

  int combinations( int n, 
                    int k );

  int factorial( int n );

protected:
  virtual void postProcessInput() override final;

  virtual void postRestartInitialization() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  std::vector< array2d< localIndex > > m_mappedNodes; // mappedNodes[subregion index][particle index][node index]. dims = {# of subregions,
                                                      // # of particles, # of nodes a particle on the subregion maps to}
  std::vector< array2d< real64 > > m_shapeFunctionValues; // mappedNodes[subregion][particle][nodal shape function value]. dims = {# of
                                                          // subregions, # of particles, # of nodes a particle on the subregion maps to}
  std::vector< array3d< real64 > > m_shapeFunctionGradientValues; // mappedNodes[subregion][particle][nodal shape function gradient
                                                                  // value][direction]. dims = {# of subregions, # of particles, # of nodes
                                                                  // a particle on the subregion maps to, 3}

  int m_solverProfiling;
  std::vector< real64 > m_profilingTimes;
  std::vector< std::string > m_profilingLabels;

  array1d< string > m_plottableFields;
  SortedArray< string > m_plottableFieldsSorted;

  TimeIntegrationOption m_timeIntegrationOption;
  UpdateMethodOption m_updateMethod;
  int m_updateOrder;

  MPI_iCommData m_iComm;

  int m_prescribedBcTable;
  array1d< int > m_boundaryConditionTypes; // TODO: Surely there's a way to have just one variable here
  array2d< real64 > m_bcTable;

  int m_prescribedFTable;
  int m_prescribedBoundaryFTable;
  InterpolationOption m_fTableInterpType;
  array2d< real64 > m_fTable;
  array1d< real64 > m_domainF;
  array1d< real64 > m_domainL;

  int m_enablePrescribedBoundaryTransverseVelocities;
  array2d< real64 > m_prescribedBoundaryTransverseVelocities; // 2 in-plane directions * 6 faces 

  array1d< real64 > m_globalFaceReactions;

  array1d< real64 > m_bodyForce;

  array1d< int > m_stressControl;
  InterpolationOption m_stressTableInterpType;
  array2d< real64 > m_stressTable;
  real64 m_stressControlKp;
  real64 m_stressControlKi;
  real64 m_stressControlKd;
  array1d< real64 > m_domainStress;
  array1d< real64 > m_stressControlLastError;
  array1d< real64 > m_stressControlITerm;

  int m_boxAverageHistory;
  real64 m_boxAverageWriteInterval;
  real64 m_nextBoxAverageWriteTime;
  array1d< real64 > m_boxAverageMin;
  array1d< real64 > m_boxAverageMax;

  int m_reactionHistory;
  real64 m_reactionWriteInterval;
  real64 m_nextReactionWriteTime;

  real64 m_explicitSurfaceNormalInfluence;
  int m_computeSurfaceNormalsOnlyOnInitialization;

  // Cohesive law variables
  int m_referenceCohesiveZone;
  int m_enableCohesiveLaws;
  CohesiveLawOption m_cohesiveLaw;
  // int m_cohesiveFieldPartitioning;
  int m_enableCohesiveFailure;
  int m_preventCZInterpentration;
  real64 m_normalForceConstant;
  real64 m_shearForceConstant;

  real64 m_numSurfaceIntegrationPoints;
  real64 m_maxCohesiveNormalStress;
  real64 m_maxCohesiveShearStress;
  real64 m_characteristicNormalDisplacement;
  real64 m_characteristicTangentialDisplacement;
  real64 m_maxCohesiveNormalDisplacement;
  real64 m_maxCohesiveTangentialDisplacement;
  SortedArray< globalIndex >  m_cohesiveNodeGlobalIndices;
  array2d< real64 > m_referenceCohesiveGridNodePartitioningSurfaceNormals;
  array2d< real64 > m_referenceCohesiveGridNodeAreas;
  array2d< real64 > m_referenceCohesiveGridNodePositions;
  array2d< real64 > m_referenceCohesiveGridNodeDamages;
  array3d< real64 > m_referenceCohesiveGridNodeSurfaceNormals;
  array3d< real64 > m_maxCohesiveGridNodeNormalDisplacement;
  array3d< real64 > m_maxCohesiveGridNodeTangentialDisplacement;

  int m_needsNeighborList;
  real64 m_neighborRadius;
  int m_binSizeMultiplier;
  real64 m_thinFeatureDFGThreshold;

  int m_useDamageAsSurfaceFlag;

  int m_FSubcycles;
  int m_LBar;
  real64 m_LBarScale;
  int m_exactJIntegration;
  real64 m_maxParticleVelocity;
  real64 m_maxParticleVelocitySquared;
  real64 m_minParticleJacobian;
  real64 m_maxParticleJacobian;
  
  OverlapCorrectionOption m_overlapCorrection;
  real64 m_overlapThreshold1;
  real64 m_overlapThreshold2;
  int m_computeSPHJacobian;

  // Currently initializes all particles to this temperature
  // TODO: read in from particle file
  real64 m_initialTemperature;
  int m_shockHeating;
  int m_useArtificialViscosity;
  real64 m_artificialViscosityQ0;
  real64 m_artificialViscosityQ1;

  int m_cpdiDomainScaling;

  real64 m_smallMass;

  int m_numContactGroups, m_numContactFlags, m_numVelocityFields;
  real64 m_separabilityMinDamage;
  int m_treatFullyDamagedAsSingleField;
  int m_surfaceDetection;
  int m_damageFieldPartitioning;
  int m_useSurfacePositionForContact;
  ContactNormalTypeOption m_contactNormalType;
  ContactGapCorrectionOption m_contactGapCorrection;
  // int m_directionalOverlapCorrection;

  int m_resetDefGradForFullyDamagedParticles;
  int m_plotUnscaledParticles;

  real64 m_frictionCoefficient;
  array2d< real64 > m_frictionCoefficientTable; 

  int m_planeStrain;
  int m_numDims;

  int m_generalizedVortexMMS;

  array1d< real64 > m_hEl;                // Grid spacing in x-y-z
  array1d< real64 > m_xLocalMin;          // Minimum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMax;          // Maximum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMinNoGhost;   // Minimum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xLocalMaxNoGhost;   // Maximum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xGlobalMin;         // Minimum global grid coordinate excluding buffer nodes
  array1d< real64 > m_xGlobalMax;         // Maximum global grid coordinate excluding buffer nodes
  array1d< real64 > m_partitionExtent;    // Length of each edge of partition including buffer and ghost cells
  array1d< real64 > m_domainExtent;       // Length of each edge of global domain excluding buffer cells
  array1d< int > m_nEl;                   // Number of elements in each grid direction including buffer and ghost cells

  // //Initial grid parameters
  // array1d< real64 > m_hEl0;
  // array1d< real64 > m_xLocalMin0;
  // array1d< real64 > m_xLocalMax0;

  array3d< int > m_ijkMap;        // Map from indices in each spatial dimension to local node ID

  int m_useEvents;                   // Events flag
  MPMEventManager* m_mpmEventManager;

  int m_surfaceHealing;

  int m_debugFlag;

  int m_computeXProfile;
  real64 m_xProfileWriteInterval;
  real64 m_nextXProfileWriteTime;
  real64 m_xProfileVx0;

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

  virtual void setConstitutiveNames( ParticleSubRegionBase & subRegion ) const override;
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
              "Symmetry" );

ENUM_STRINGS( SolidMechanicsMPM::InterpolationOption,
              "Linear",
              "Cosine",
              "Smoothstep" );

ENUM_STRINGS( SolidMechanicsMPM::ContactNormalTypeOption,
              "Difference",
              "MassWeighted",
              "LargerMass",
              "Mixed" );

ENUM_STRINGS( SolidMechanicsMPM::ContactGapCorrectionOption,
              "Simple",
              "Implicit",
              "Softened" );

ENUM_STRINGS( SolidMechanicsMPM::OverlapCorrectionOption,
              "Off",
              "NormalForce",
              "SPH" );

ENUM_STRINGS( SolidMechanicsMPM::CohesiveLawOption,
              "Uncoupled",
              "PPR",
              "NeedlemanXu" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


} /* namespace geos */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
