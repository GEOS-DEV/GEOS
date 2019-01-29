/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FlowSolverBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_FLOWSOLVERBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_FLOWSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class FlowSolverBase
 *
 * Base class for finite volume fluid flow solvers.
 * Provides some common features
 */
class FlowSolverBase : public SolverBase
{
public:
/**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  FlowSolverBase( const std::string& name,
                  ManagedGroup * const parent );


  /// deleted default constructor
  FlowSolverBase() = delete;

  /// deleted copy constructor
  FlowSolverBase( FlowSolverBase const & ) = delete;

  /// default move constructor
  FlowSolverBase( FlowSolverBase && ) = default;

  /// deleted assignment operator
  FlowSolverBase & operator=( FlowSolverBase const & ) = delete;

  /// deleted move operator
  FlowSolverBase & operator=( FlowSolverBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~FlowSolverBase() override;

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBodies ) override;

  void setPoroElasticCoupling() { m_poroElasticFlag = 1; }

  localIndex fluidIndex() const { return m_fluidIndex; }

  localIndex solidIndex() const { return m_solidIndex; }

  localIndex numDofPerCell() const { return m_numDofPerCell; }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // input data
    static constexpr auto referencePorosityString = "referencePorosity";
    static constexpr auto permeabilityString      = "permeability";

    // gravity term precomputed values
    static constexpr auto gravityFlagString  = "gravityFlag";
    static constexpr auto gravityDepthString = "gravityDepth";

    // misc inputs
    static constexpr auto fluidNameString      = "fluidName";
    static constexpr auto solidNameString      = "solidName";
    static constexpr auto fluidIndexString     = "fluidIndex";
    static constexpr auto solidIndexString     = "solidIndex";

    using ViewKey = dataRepository::ViewKey;

    // input data
    ViewKey referencePorosity = { referencePorosityString };
    ViewKey permeability      = { permeabilityString };

    // gravity term precomputed values
    ViewKey gravityFlag  = { gravityFlagString };
    ViewKey gravityDepth = { gravityDepthString };

    // misc inputs
    ViewKey discretization = { discretizationString };
    ViewKey fluidName      = { fluidNameString };
    ViewKey solidName      = { solidNameString };
    ViewKey fluidIndex     = { fluidIndexString };
    ViewKey solidIndex     = { solidIndexString };

  } viewKeysFlowSolverBase;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysFlowSolverBase;

private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void PrecomputeData(DomainPartition *const domain);

protected:

  virtual void InitializePreSubGroups(ManagedGroup * const rootGroup) override;

  virtual void InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup) override;


  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void ResetViews( DomainPartition * const domain );

  /// flag to determine whether or not to apply gravity
  integer m_gravityFlag;

  /// name of the fluid constitutive model
  string m_fluidName;

  /// name of the solid constitutive model
  string m_solidName;

  /// index of the fluid constitutive model
  localIndex m_fluidIndex;

  /// index of the solid constitutive model
  localIndex m_solidIndex;

  /// flag to determine whether or not coupled with solid solver
  integer m_poroElasticFlag;

  /// the number of Degrees of Freedom per cell
  localIndex m_numDofPerCell;

  /// views into constant data fields
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> m_elemGhostRank;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  m_volume;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  m_gravDepth;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  m_porosityRef;

};

}

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_FLOWSOLVERBASE_HPP_
