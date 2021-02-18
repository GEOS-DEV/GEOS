/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_

#include "codingUtilities/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers
#include "managers/FieldSpecification/FieldSpecificationManager.hpp" // a manager that can access and set values on the discretized domain
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"  // interface to linear solvers and linear algebra libraries

namespace geosx
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

//START_SPHINX_INCLUDE_02
class LaplaceFEM : public SolverBase
{
public:
  // The default nullary constructor is disabled to avoid compiler auto-generation:
  LaplaceFEM() = delete;

  // The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  LaplaceFEM( const string & name,
              Group * const parent );

  // Destructor
  virtual ~LaplaceFEM() override;

  // "CatalogName()" return the string used as XML tag in the input file.
  // It ties the XML tag with this C++ classes. This is important.
  static string catalogName() { return "LaplaceFEM"; }

  // This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override final;

//END_SPHINX_INCLUDE_02
/**
 * @defgroup Solver Interface Functions
 *
 * These functions provide the primary interface that is required for derived classes
 */
/**@{*/

  //START_SPHINX_INCLUDE_03
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = false ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
    resetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  //END_SPHINX_INCLUDE_03
  /**@}*/

  // This method is specific to this Laplace solver
  // It is used to apply Dirichlet boundary condition
  // and called when the base class applyBoundaryConditions() is called
  void applyDirichletBCImplicit( real64 const time,
                                 DofManager const & dofManager,
                                 DomainPartition & domain,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

  // Choice of transient treatment options (steady, backward, forward Euler scheme):
  //START_SPHINX_INCLUDE_01
  enum class TimeIntegrationOption : integer
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };
  //END_SPHINX_INCLUDE_01


  // This structure stores ``dataRepository::ViewKey`` objects
  // used as binding between the input XML tags and source
  // code variables (here, timeIntegrationOption and fieldVarName)
  //START_SPHINX_INCLUDE_04
  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
  } laplaceFEMViewKeys;
  //END_SPHINX_INCLUDE_04

private:

  // These two classes are specific to the Laplace solver:
  string m_fieldName;  // User-defined name of the physical quantity we wish to solve for (such as "Temperature", etc.)
  TimeIntegrationOption m_timeIntegrationOption;  // Choice of transient treatment (SteadyState, ImplicitTransient or ExplicitTransient)

};


/* REGISTERING NEW ENUMS:
   --------------------------
   In order to register an enumeration type with the Data Repository and have its value read from input,
   we must define stream insertion/extraction operators. This is a common task, so GEOSX provides
   a facility for automating it. Upon including ``common/EnumStrings.hpp``, we can call the following macro
   at the namespace scope (in this case, right after the ``LaplaceFEM`` class definition is complete):
 */
//START_SPHINX_INCLUDE_05
ENUM_STRINGS( LaplaceFEM::TimeIntegrationOption, "SteadyState", "ImplicitTransient", "ExplicitTransient" )
//END_SPHINX_INCLUDE_05

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_ */
