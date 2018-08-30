/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
class BoundaryConditionBase;
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

  virtual void FillDocumentationNode() override;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup ) override;

  virtual void FinalInitialization( dataRepository::ManagedGroup * const rootGroup ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {

    // input data
    constexpr static auto referencePorosityString = "referencePorosity";
    constexpr static auto permeabilityString      = "permeability";

    // gravity term precomputed values
    constexpr static auto gravityFlagString  = "gravityFlag";
    constexpr static auto gravityDepthString = "gravityDepth";

    // misc inputs
    constexpr static auto discretizationString = "discretization";
    constexpr static auto fluidNameString      = "fluidName";
    constexpr static auto solidNameString      = "solidName";
    constexpr static auto fluidIndexString     = "fluidIndex";
    constexpr static auto solidIndexString     = "solidIndex";

  } viewKeys;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeys;

private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void PrecomputeData(DomainPartition *const domain);

protected:

  /// flag to determine whether or not to apply gravity
  integer m_gravityFlag;

  /// name of the FV discretization object in the data repository
  string m_discretizationName;

  /// name of the fluid constitutive model
  string m_fluidName;

  /// name of the solid constitutive model
  string m_solidName;

  /// index of the fluid constitutive model
  localIndex m_fluidIndex;

  /// index of the solid constitutive model
  localIndex m_solidIndex;

};

}

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_FLOWSOLVERBASE_HPP_
