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
 * @file MatrixFreeSolidMechanicsFEMOperator.hpp
 */


#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MATRICFREESOLIDMECHANICSFEMOPERATOR_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MATRICFREESOLIDMECHANICSFEMOPERATOR_HPP_


#include "dataRepository/ExecutableGroup.hpp"
#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/common/LinearOperatorWithBC.hpp"

namespace geos
{

class MatrixFreeSolidMechanicsFEMOperator : public LinearOperator< ParallelVector >
{
private:
  dataRepository::Group & m_meshBodies;
  map< std::pair< string, string >, array1d< string > > const & m_meshTargets;
  DofManager & m_dofManager;
  string const m_finiteElementName;
  int const m_kernelOptimizationOption = 0;

public:
  MatrixFreeSolidMechanicsFEMOperator( DomainPartition & domain, 
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets, 
                                       DofManager & dofManager, 
                                       string const & finiteElementName,
                                       int const kernelOptimizationOption );

  MatrixFreeSolidMechanicsFEMOperator( dataRepository::Group & meshBodies, 
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets, 
                                       DofManager & dofManager,
                                       string const & finiteElementName,
                                       int const kernelOptimizationOption );

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const;

  void computeDiagonal( ParallelVector & diagonal ) const;

  virtual globalIndex numGlobalRows() const;

  virtual globalIndex numGlobalCols() const;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  virtual MPI_Comm comm() const;
};

} /* namespace geos */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MATRICFREESOLIDMECHANICSFEMOPERATOR_HPP_ */
