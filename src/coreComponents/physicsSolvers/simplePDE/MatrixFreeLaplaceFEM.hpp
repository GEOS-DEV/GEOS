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

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_EXPLICIT_LAPLACE_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_EXPLICIT_LAPLACE_FEM_HPP_

#include "dataRepository/ExecutableGroup.hpp"
#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"

namespace geosx
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

class MatrixFreeLaplaceFEMOperator : public LinearOperator< ParallelVector >
{
private:
  dataRepository::Group & m_meshBodies;
  map< string, array1d< string > > m_meshTargets;
  DofManager & m_dofManager;

public:
  MatrixFreeLaplaceFEMOperator( DomainPartition & domain, DofManager & dofManager );
  MatrixFreeLaplaceFEMOperator( dataRepository::Group & meshBodies, DofManager & dofManager );

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const;

  virtual globalIndex numGlobalRows() const;

  virtual globalIndex numGlobalCols() const;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  virtual MPI_Comm comm() const;
};

template < typename Vector >
class LinearOperatorWithBC : public LinearOperator< Vector >
{
private:
  LinearOperator< Vector > const & m_unconstrained_op;

public:
  LinearOperatorWithBC( LinearOperator< Vector > const & unconstrained_op, DomainPartition & domain );

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const;

  virtual globalIndex numGlobalRows() const;

  virtual globalIndex numGlobalCols() const;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  virtual MPI_Comm comm() const;
};

class MatrixFreePreconditionerIdentity : public PreconditionerIdentity< HypreInterface >
{
private:
  DofManager & m_dofManager;

public:
  MatrixFreePreconditionerIdentity( DofManager & dofManager );

  virtual globalIndex numGlobalRows() const;

  virtual globalIndex numGlobalCols() const;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  virtual MPI_Comm comm() const;
};

//START_SPHINX_INCLUDE_BEGINCLASS
class MatrixFreeLaplaceFEM : public LaplaceBaseH1
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  MatrixFreeLaplaceFEM() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  MatrixFreeLaplaceFEM( const string & name,
                      Group * const parent );

  /// Destructor
  virtual ~MatrixFreeLaplaceFEM() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "MatrixFreeLaplaceFEM"; }

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

};
} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_EXPLICIT_LAPLACE_FEM_HPP_ */
