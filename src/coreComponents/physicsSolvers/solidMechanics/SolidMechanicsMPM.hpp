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
    static constexpr char const * vTildeString() { return "velocityTilde"; }
    static constexpr char const * uhatTildeString() { return "uhatTilde"; }
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * newmarkGammaString() { return "newmarkGamma"; }
    static constexpr char const * newmarkBetaString() { return "newmarkBeta"; }
    static constexpr char const * massDampingString() { return "massDamping"; }
    static constexpr char const * stiffnessDampingString() { return "stiffnessDamping"; }
    static constexpr char const * useVelocityEstimateForQSString() { return "useVelocityForQS"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * maxNumResolvesString() { return "maxNumResolves"; }
    static constexpr char const * strainTheoryString() { return "strainTheory"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * forceExternalString() { return "externalForce"; }
    static constexpr char const * contactRelationNameString() { return "contactRelationName"; }
    static constexpr char const * noContactRelationNameString() { return "NOCONTACT"; }
    static constexpr char const * contactForceString() { return "contactForce"; }
    static constexpr char const * maxForceString() { return "maxForce"; }
    static constexpr char const * elemsAttachedToSendOrReceiveNodesString() { return "elemsAttachedToSendOrReceiveNodes"; }
    static constexpr char const * elemsNotAttachedToSendOrReceiveNodesString() { return "elemsNotAttachedToSendOrReceiveNodes"; }

    static constexpr char const * sendOrReceiveNodesString() { return "sendOrReceiveNodes";}
    static constexpr char const * nonSendOrReceiveNodesString() { return "nonSendOrReceiveNodes";}
    static constexpr char const * targetNodesString() { return "targetNodes";}


    dataRepository::ViewKey vTilde = { vTildeString() };
    dataRepository::ViewKey uhatTilde = { uhatTildeString() };
    dataRepository::ViewKey newmarkGamma = { newmarkGammaString() };
    dataRepository::ViewKey newmarkBeta = { newmarkBetaString() };
    dataRepository::ViewKey massDamping = { massDampingString() };
    dataRepository::ViewKey stiffnessDamping = { stiffnessDampingString() };
    dataRepository::ViewKey useVelocityEstimateForQS = { useVelocityEstimateForQSString() };
    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;


  SortedArray< localIndex > & getElemsAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() );
  }

  SortedArray< localIndex > & getElemsNotAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() );
  }

  real64 & getMaxForce() { return m_maxForce; }

  arrayView1d< ParallelVector > const & getRigidBodyModes() const
  {
    return m_rigidBodyModes;
  }

  array1d< ParallelVector > & getRigidBodyModes()
  {
    return m_rigidBodyModes;
  }

  std::vector<int> getNodes(std::vector<int> const & cellID);

  std::vector<real64> getWeights(LvArray::ArraySlice<double, 1, 0, long> const & p_x,
                                 std::vector<int> const & cellID,
                                 arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & g_X);

  void initialize(arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X,
                  ParticleManager & particleManager);

protected:
  virtual void postProcessInput() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

//  virtual void initializePostInitialConditionsPreSubGroups() override final;

  real64 m_newmarkGamma;
  real64 m_newmarkBeta;
  real64 m_massDamping;
  real64 m_stiffnessDamping;
  TimeIntegrationOption m_timeIntegrationOption;
  integer m_useVelocityEstimateForQS;
  real64 m_maxForce = 0.0;
  integer m_maxNumResolves;
  integer m_strainTheory = 0;
  string m_contactRelationName;
  MPI_iCommData m_iComm;

  std::vector<real64> m_hx{DBL_MAX,DBL_MAX,DBL_MAX}; // Grid spacing in x-y-z
  std::vector<real64> m_xMin{DBL_MAX,DBL_MAX,DBL_MAX}; // Minimum grid coordinate
  std::vector<real64> m_xMax{DBL_MIN,DBL_MIN,DBL_MIN}; // Maximum grid coordinate
  std::vector<real64> m_domainL{0.0,0.0,0.0}; // Length of each edge of grid
  std::vector<int> m_nEl{0,0,0}; // Number of elements in each grid direction
  std::vector<std::vector<std::vector<int>>> m_ijkMap;

  /// Rigid body modes
  array1d< ParallelVector > m_rigidBodyModes;

private:
  virtual void setConstitutiveNames( ParticleSubRegionBase & subRegion ) const override;

};

ENUM_STRINGS( SolidMechanicsMPM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
