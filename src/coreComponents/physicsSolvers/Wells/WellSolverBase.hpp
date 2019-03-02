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
 * @file WellSolverBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_WELLS_WELLSOLVERBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_WELLS_WELLSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}
class DomainPartition;

/**
 * @class WellSolverBase
 *
 * Base class for well solvers.
 * Provides some common features
 */
class WellSolverBase : public SolverBase
{
public:

  struct SubRegionTag
  {
    static constexpr integer RES  = 0;
    static constexpr integer WELL = 1;
  };

  struct ElemTag
  {
    static constexpr integer PREV = 0;
    static constexpr integer NEXT = 1;
  };

/**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  WellSolverBase( const std::string& name,
                  ManagedGroup * const parent );

  /// deleted default constructor
  WellSolverBase() = delete;

  /// deleted copy constructor
  WellSolverBase( WellSolverBase const & ) = delete;

  /// default move constructor
  WellSolverBase( WellSolverBase && ) = default;

  /// deleted assignment operator
  WellSolverBase & operator=( WellSolverBase const & ) = delete;

  /// deleted move operator
  WellSolverBase & operator=( WellSolverBase && ) = delete;

  localIndex fluidIndex() const { return m_fluidIndex; }

  localIndex numDofPerElement() const { return m_numDofPerElement; }

  localIndex numDofPerConnection() const { return m_numDofPerConnection; }

  localIndex numDofPerResElement() const { return m_numDofPerResElement; }

  globalIndex getFirstWellElementDofNumber() const { return m_firstWellElemDofNumber; }

  globalIndex getElementOffset( globalIndex welemDofNumber ) const;
  
  /**
   * @brief default destructor
   */
  virtual ~WellSolverBase() override;

  virtual void RegisterDataOnMesh(ManagedGroup * const meshBodies) override;

  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
  virtual void SetSparsityPattern( DomainPartition const * const domain,
                                   Epetra_FECrsGraph * const sparsity,
				   globalIndex firstWellElemDofNumber,
				   localIndex numDofPerResElement );

  /**
   * @brief sets the dof indices for this solver
   * @param meshLevel the mesh object (single level only)
   * @param numLocalRows the number of local rows on this partition
   * @param numGlobalRows the number of global rows in the problem
   * @param offset the DOF offset for this solver in the case of a non-block system
   *
   * This function sets the number of global rows, and sets the dof numbers for
   * this solver. dof numbers are referred to trilinosIndices currently.
   */
  virtual void SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                             localIndex & numLocalRows,
                                             globalIndex & numGlobalRows,
                                             localIndex offset );

  
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // gravity term precomputed values
    static constexpr auto gravityFlagString  = "gravityFlag";
    static constexpr auto gravityDepthString = "gravityDepth";

    // misc inputs
    static constexpr auto fluidNameString  = "fluidName";
    static constexpr auto fluidIndexString = "fluidIndex";

    // bhp control
    static constexpr auto bhpString = "BHP";

    using ViewKey = dataRepository::ViewKey;

    // gravity term precomputed values
    ViewKey gravityFlag  = { gravityFlagString };
    ViewKey gravityDepth = { gravityDepthString };

    // misc inputs
    ViewKey fluidName  = { fluidNameString };
    ViewKey fluidIndex = { fluidIndexString };

    // bhp control
    ViewKey bhp = { bhpString }; 
  
  } viewKeysWellSolverBase;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysWellSolverBase;

private:
  
  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void PrecomputeData(DomainPartition *const domain);
  
protected:

  virtual void InitializePreSubGroups(ManagedGroup * const rootGroup) override;

  virtual void InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup) override;

  virtual void ResetViews( DomainPartition * const domain );
  
  /// flag to determine whether or not to apply gravity
  integer m_gravityFlag;

  /// name of the fluid constitutive model
  string m_fluidName;

  /// index of the fluid constitutive model
  localIndex m_fluidIndex;

  /// the number of Degrees of Freedom per well element
  localIndex m_numDofPerElement;

  /// the number of Degrees of Freedom per well connection
  localIndex m_numDofPerConnection;

  /// the number of Degrees of Freedom per reservoir element
  localIndex m_numDofPerResElement;

  /// the number of the first Degree of Freedom corresponding to a well var/eq
  globalIndex m_firstWellElemDofNumber;

  /// flag to determine whether the well mass balance equations are normalized
  integer m_normalizeMassBalanceEqnsFlag;
  
  /// views into reservoir constant data fields
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  m_resGravDepth;
};

}

#endif //SRC_COMPONENTS_CORE_SRC_WELLS_WELLSOLVERBASE_HPP_
