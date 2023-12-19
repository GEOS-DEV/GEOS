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

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_BASE_HPP
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_BASE_HPP

#include "MPMSolverFields.hpp"
// #include "codingUtilities/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers
// #include "fieldSpecification/FieldSpecificationManager.hpp" // a manager that can access and set values on the discretized domain

namespace geos
{
    
class MPMBase : public SolverBase
{
public:
  /**
   * Nullary Constructor
   * The default nullary constructor is disabled to avoid compiler auto-generation:
  */
  MPMBase() = delete;

  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  MPMBase( const string & name,
                 Group * const parent );             

  /**
   * Destructor
  */
  virtual ~MPMBase() override;

  MPMBase( MPMBase const & ) = delete;
  MPMBase( MPMBase && ) = default;    

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * massString() { return "mass"; }
    static constexpr char const * velocityString() { return "velocity"; }
    static constexpr char const * momentumString() { return "momentum"; }
    static constexpr char const * accelerationString() { return "acceleration"; }
    static constexpr char const * forceExternalString() { return "externalForce"; }
    static constexpr char const * forceInternalString() { return "internalForce"; }
    static constexpr char const * forceContactString() { return "contactForce"; }
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * damageGradientString() { return "damageGradient"; }
    static constexpr char const * maxDamageString() { return "maxDamage"; }
    static constexpr char const * surfaceNormalString() { return "surfaceNormal"; }
    static constexpr char const * materialPositionString() { return "materialPosition"; }
    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }
  } solidMechanicsViewKeys;

  virtual void initializePreSubGroups() override;

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  //START_SPHINX_INCLUDE_SOLVERINTERFACE
  // virtual real64 solverStep( real64 const & time_n, real64 const & dt, integer const cycleNumber, DomainPartition & domain ) override;

  // virtual real64 explicitStep( real64 const & time_n, real64 const & dt, integer const cycleNumber, DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override; // TODO: maybe final;
  //END_SPHINX_INCLUDE_SOLVERINTERFACE

  // // This method is called when its host event is triggered
  // virtual bool execute( real64 const time_n, real64 const dt, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;

  //START_SPHINX_INCLUDE_REGENUM
  /**
   * @enum TimeIntegrationOption
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
   * The options for essential boundary conditions
   */
  enum struct BoundaryConditionOption : integer
  {
    OUTFLOW,    //!<Outflow
    SYMMETRY    //!<Symmetry
  };
  //END_SPHINX_INCLUDE_REGENUM

protected:

  int m_solverProfiling;
  std::vector< real64 > m_profilingTimes;
  std::vector< std::string > m_profilingLabels;
  int m_damageFieldPartitioning;
  int m_surfaceDetection;
  int m_directionalOverlapCorrection; 
  int m_needsNeighborList;
  int m_numDims;
  int m_planeStrain;
  array1d< int > m_boundaryConditionTypes; // TODO: Surely there's a way to have just one variable here
 

  virtual void postProcessInput() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  virtual void setConstitutiveNames( ParticleSubRegionBase & subRegion ) const override;

private:

};

//START_SPHINX_INCLUDE_REGENUM
///* REGISTERING NEW ENUMS:
ENUM_STRINGS( MPMBase::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

ENUM_STRINGS( MPMBase::BoundaryConditionOption,
              "Outflow",
              "Symmetry" );
//END_SPHINX_INCLUDE_REGENUM

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_BASE_HPP