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
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "physicsSolvers/SolverBase.hpp"

#include "SolidMechanicsLagrangianFEMKernels.hpp"

namespace geosx
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
  virtual ~SolidMechanicsMPM() override;

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
    GEOSX_UNUSED_VAR( domain );
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
    static constexpr char const * momentumString() { return "momentum"; }
    static constexpr char const * forceContactString() { return "contactForce"; }
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * maxDamageString() { return "maxDamage"; }
    static constexpr char const * surfaceNormalString() { return "surfaceNormal"; }
    static constexpr char const * materialPositionString() { return "materialPosition"; }

    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }

    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  void initialize (NodeManager & nodeManager,
                   ParticleManager & particleManager,
                   SpatialPartition & partition );

protected:
  virtual void postProcessInput() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  int m_solverProfiling;

  std::vector< real64 > m_profilingTimes;
  std::vector< std::string > m_profilingLabels;

  TimeIntegrationOption m_timeIntegrationOption;
  MPI_iCommData m_iComm;

  array1d< int > m_boundaryConditionTypes;

  real64 m_smallMass;

  int m_numContactGroups, m_numContactFlags, m_numVelocityFields; // TODO: I think we only need m_numVelocityFields
  int m_damageFieldPartitioning;
  int m_contactGapCorrection;
  real64 m_frictionCoefficient;

  int m_planeStrain;
  int m_numDims;



  std::array< real64, 3 > m_hEl;                // Grid spacing in x-y-z
  std::array< real64, 3 > m_xLocalMin;          // Minimum local grid coordinate including ghost nodes
  std::array< real64, 3 > m_xLocalMax;          // Maximum local grid coordinate including ghost nodes
  std::array< real64, 3 > m_xLocalMinNoGhost;   // Minimum local grid coordinate EXCLUDING ghost nodes
  std::array< real64, 3 > m_xLocalMaxNoGhost;   // Maximum local grid coordinate EXCLUDING ghost nodes
  std::array< real64, 3 > m_xGlobalMin;         // Minimum global grid coordinate excluding buffer nodes
  std::array< real64, 3 > m_xGlobalMax;         // Maximum global grid coordinate excluding buffer nodes
  std::array< real64, 3 > m_domainLengths;      // Length of each edge of grid
  std::array< int, 3 > m_nEl;                   // Number of elements in each grid direction including buffer cells
  array3d< int > m_ijkMap;                      // Map from cell-spaced coordinates to cell ID

  int m_voigtMap[3][3];

private:
  virtual void setConstitutiveNames( ParticleSubRegionBase & subRegion ) const override;

  void syncGridFields( std::vector< std::string > const & fieldNames,
                       DomainPartition & domain,
                       NodeManager & nodeManager,
                       MeshLevel & mesh );

  void enforceGridVectorFieldSymmetryBC( array3d< real64 > & vectorField,
                                         arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridReferencePosition,
                                         Group & nodeSets );

  void computeGridSurfaceNormals( ParticleManager & particleManager,
                                  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridReferencePosition,
                                  array3d< real64 > & gridSurfaceNormal );

  void normalizeGridSurfaceNormals( array2d< real64 > & gridMass,
                                    array3d< real64 > & gridSurfaceNormal );

  void computeContactForces( const real64 dt,
                             const array2d< real64 > & gridMass,
                             const array2d< real64 > & GEOSX_UNUSED_PARAM( gridDamage ),
                             const array2d< real64 > & GEOSX_UNUSED_PARAM( gridMaxDamage ),
                             const array3d< real64 > & gridVelocity,
                             const array3d< real64 > & gridMomentum,
                             const array3d< real64 > & gridSurfaceNormal,
                             const array3d< real64 > & gridMaterialPosition,
                             array3d< real64 > & gridContactForce );

  void computePairwiseNodalContactForce( const int & separable,
                                         const real64 & dt,
                                         const real64 & mA,
                                         const real64 & mB,
                                         const arraySlice1d< real64 > vA,
                                         const arraySlice1d< real64 > GEOSX_UNUSED_PARAM( vB ),
                                         const arraySlice1d< real64 > qA,
                                         const arraySlice1d< real64 > qB,
                                         const arraySlice1d< real64 > nA,
                                         const arraySlice1d< real64 > nB,
                                         const arraySlice1d< real64 > xA, // Position of field A
                                         const arraySlice1d< real64 > xB, // Position of field B
                                         arraySlice1d< real64 > fA,
                                         arraySlice1d< real64 > fB );

  real64 subtractDot( const arraySlice1d< real64 > & u,
                      const arraySlice1d< real64 > & v,
                      const arraySlice1d< real64 > & n );

  void computeOrthonormalBasis( const array1d< real64 > & e1,  // input "normal" unit vector.
                                array1d< real64 > & e2,        // output "tangential" unit vector.
                                array1d< real64 > & e3 );      // output "tangential" unit vector.

  void setGridFieldLabels( NodeManager & nodeManager );

  void solverProfiling( std::string label );

};

ENUM_STRINGS( SolidMechanicsMPM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

ENUM_STRINGS( SolidMechanicsMPM::BoundaryConditionOption,
              "Outflow",
              "Symmetry" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
