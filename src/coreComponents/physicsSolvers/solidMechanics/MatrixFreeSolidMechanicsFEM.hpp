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

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_

#include "dataRepository/ExecutableGroup.hpp"
#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"

namespace geosx
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

class MatrixFreeSolidMechanicsFEMOperator : public LinearOperator< ParallelVector >
{
private:
  dataRepository::Group & m_meshBodies;
  map< string, array1d< string > > & m_meshTargets;
  DofManager & m_dofManager;
  string const & m_finiteElementName;

public:
  MatrixFreeSolidMechanicsFEMOperator( DomainPartition & domain, map< string, array1d< string > > & meshTargets, DofManager & dofManager, string const & finiteElementName );
  MatrixFreeSolidMechanicsFEMOperator( dataRepository::Group & meshBodies, map< string, array1d< string > > & meshTargets, DofManager & dofManager, string const & finiteElementName );

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const;

  void computeDiagonal( ParallelVector & diagonal ) const;

  virtual globalIndex numGlobalRows() const;

  virtual globalIndex numGlobalCols() const;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  virtual MPI_Comm comm() const;
};

template < typename Vector >
class LinearOperatorWithBC : public LinearOperator< Vector >
{
public:
  enum class DiagPolicy : integer
  {
    KeepDiagonal, ///< Use diagonal values
    DiagonalOne,  ///< Use one on the diagonal
    DiagonalZero, ///< Use zero on the diagonal
  };

  LinearOperatorWithBC( SolverBase const & solver,
                        LinearOperator< Vector > const & unconstrained_op,
                        DomainPartition & domain,
                        DofManager const & dofManager,
                        string fieldName,
                        real64 const time,
                        DiagPolicy diagPolicy = DiagPolicy::DiagonalOne ):
    m_unconstrained_op( unconstrained_op ),
    m_domain( domain ),
    m_dofManager( dofManager ),
    m_fieldName( fieldName ),
    m_time( time ),
    m_diagPolicy( diagPolicy ),
    m_diagonal( m_dofManager.numLocalDofs() )
  {
    using POLICY = parallelDevicePolicy<>;
    switch (m_diagPolicy)
    {
      case DiagPolicy::DiagonalOne:
        {
          m_diagonal.setValues< POLICY >( 1.0 );
        }
        break;
      case DiagPolicy::DiagonalZero:
        {
          m_diagonal.zero();
        }
        break;
      case DiagPolicy::KeepDiagonal:
        GEOSX_ERROR("Not yet implemented.");
        // TODO: compute diagonal
        break;
    }
    // compute m_constrainedDofIndices and m_rhsContributions
    FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();
    globalIndex totalSize = 0;
    solver.forMeshTargets( m_domain.getMeshBodies(), [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
    {
      fsManager.apply( m_time,
                       mesh,
                       "nodeManager",
                       m_fieldName,
                       [&]( FieldSpecificationBase const & bc,
                            string const &,
                            SortedArrayView< localIndex const > const & targetSet,
                            dataRepository::Group & targetGroup,
                            string const & GEOSX_UNUSED_PARAM( fieldName ) )
      {
        totalSize += targetSet.size();
      } );
    } );
    // NOTE: we're not checking for duplicates.
    m_constrainedDofIndices.reserve( totalSize );
    m_rhsContributions.reserve( totalSize );
    solver.forMeshTargets( m_domain.getMeshBodies(), [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
    {
      auto const & nodeManager = mesh.getNodeManager();
      auto const & field = nodeManager.getReference< array1d< real64 > >( m_fieldName ).toViewConst();

      fsManager.apply( m_time,
                       mesh,
                       "nodeManager",
                       m_fieldName,
                       [&]( FieldSpecificationBase const & bc,
                            string const &,
                            SortedArrayView< localIndex const > const & targetSet,
                            dataRepository::Group & targetGroup,
                            string const & fieldName )
      {
        array1d< globalIndex > dofArray( targetSet.size() );
        arrayView1d< globalIndex > const & dof = dofArray.toView();

        array1d< real64 > rhsContributionArray( targetSet.size() );
        arrayView1d< real64 > const & rhsContribution = rhsContributionArray.toView();
        arrayView1d< globalIndex const > const & dofMap =
          targetGroup.getReference< array1d< globalIndex > >( m_dofManager.getKey( fieldName ) );
        bc.computeRhsContribution< FieldSpecificationEqual, POLICY >( targetSet,
                                                                      time,
                                                                      1.0, // TODO: double check
                                                                      targetGroup,
                                                                      dofMap,
                                                                      m_dofManager.rankOffset(),
                                                                      m_diagonal.toViewConst(),
                                                                      dof,
                                                                      rhsContribution,
                                                                      field );

        dofArray.move( LvArray::MemorySpace::host, false );
        rhsContribution.move( LvArray::MemorySpace::host, false );
        m_constrainedDofIndices.insert( m_constrainedDofIndices.size(), dofArray.begin(), dofArray.end() );
        m_rhsContributions.insert( m_rhsContributions.size(), rhsContribution.begin(), rhsContribution.end() );
      } );
    } );
    // TODO: sort m_constrainedDofIndices and m_rhsContributions together

    srcWithBC.create( dofManager.numLocalDofs(), this->comm() );
    tmpRhs.create( dofManager.numLocalDofs(), this->comm() );
  }

  void computeConstrainedRHS( ParallelVector & rhs ) const
  {
    using POLICY = parallelDevicePolicy<>;

    // Construct [x_BC,0]
    srcWithBC.zero();
    arrayView1d< real64 > const localBC = srcWithBC.open();
    arrayView1d< localIndex const > const localBCIndices = m_constrainedDofIndices.toViewConst();
    arrayView1d< real64 const > const localDiag = m_diagonal.toViewConst();
    arrayView1d< real64 const > const localRhsContributions = m_rhsContributions.toViewConst();
    forAll< POLICY >( m_constrainedDofIndices.size(),
                      [ localBC, localBCIndices, localDiag, localRhsContributions ] GEOSX_HOST_DEVICE
                        ( localIndex const i )
    {
      localIndex const idx = localBCIndices[ i ];
      localBC[ idx ] = localRhsContributions [ i ] / localDiag[ idx ];
    } );
    srcWithBC.close();

    // Bottom contribution to rhs
    tmpRhs = rhs;
    m_unconstrained_op.apply( srcWithBC, rhs );
    rhs.axpby( 1.0, tmpRhs, -1.0 );

    // D_GG x_BC
    arrayView1d< real64 > const localRhs = rhs.open();
    forAll< POLICY >( m_constrainedDofIndices.size(),
                      [ localRhs, localBCIndices, localRhsContributions ] GEOSX_HOST_DEVICE
                        ( localIndex const i )
    {
      localIndex const idx = localBCIndices[ i ];
      localRhs[ idx ] = localRhsContributions [ i ];
    } );
    rhs.close();
  }

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const
  {
    using POLICY = parallelDevicePolicy<>;

    srcWithBC = src;

    arrayView1d< real64 > const localSrcWithBC = srcWithBC.open();
    arrayView1d< localIndex const > const localBCIndices = m_constrainedDofIndices.toViewConst();
    forAll< POLICY >( m_constrainedDofIndices.size(), [localSrcWithBC, localBCIndices] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localSrcWithBC[ localBCIndices[ i ] ] = 0.0;
    } );
    srcWithBC.close();

    m_unconstrained_op.apply( srcWithBC, dst );

    arrayView1d< real64 const > const localSrc = src.values();
    arrayView1d< real64 > const localDst = dst.open();
    switch (m_diagPolicy)
    {
      case DiagPolicy::DiagonalOne:
        {
          forAll< POLICY >( m_constrainedDofIndices.size(), [ localSrc, localDst, localBCIndices ] GEOSX_HOST_DEVICE ( localIndex const i )
          {
            localIndex const idx = localBCIndices[ i ];
            localDst[ idx ] = localSrc [ idx ];
          } );
        }
        break;
      case DiagPolicy::DiagonalZero:
        {
          forAll< POLICY >( m_constrainedDofIndices.size(), [ localDst, localBCIndices ] GEOSX_HOST_DEVICE ( localIndex const i )
          {
            localIndex const idx = localBCIndices[ i ];
            localDst[ idx ] = 0.0;
          } );
        }
        break;
      case DiagPolicy::KeepDiagonal:
        GEOSX_ERROR("Not yet implemented.");
        // TODO: set dst to diag * xBC
        break;
    }
    dst.close();
  }

  globalIndex numGlobalRows() const
  {
    return m_dofManager.numGlobalDofs();
  }

  globalIndex numGlobalCols() const
  {
    return m_dofManager.numGlobalDofs();
  }

  localIndex numLocalRows() const
  {
    return m_dofManager.numLocalDofs();
  }

  localIndex numLocalCols() const
  {
    return m_dofManager.numLocalDofs();
  }

  MPI_Comm comm() const
  {
    return MPI_COMM_GEOSX;
  }

private:
  LinearOperator< Vector > const & m_unconstrained_op;
  DomainPartition & m_domain;
  DofManager const & m_dofManager;
  string m_fieldName;
  real64 const m_time;
  DiagPolicy m_diagPolicy;
  array1d< localIndex > m_constrainedDofIndices;
  array1d< real64 > m_rhsContributions;

  mutable ParallelVector srcWithBC;
  mutable ParallelVector tmpRhs;
  array1d< real64 > m_diagonal;
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
class MatrixFreeSolidMechanicsFEM : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  MatrixFreeSolidMechanicsFEM() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  MatrixFreeSolidMechanicsFEM( const string & name,
                               Group * const parent );

  /// Destructor
  virtual ~MatrixFreeSolidMechanicsFEM() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "MatrixFreeSolidMechanicsFEM"; }

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

protected:
  string m_fieldName;

};
} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_MFSOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_ */
