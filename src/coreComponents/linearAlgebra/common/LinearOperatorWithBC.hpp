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

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATORWITHBC_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATORWITHBC_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"

#include "LvArray/src/output.hpp"

namespace geos
{

template< typename T >
GEOS_HOST_DEVICE
inline real64 bcFieldValue( T const & field, 
                            localIndex const index, 
                            int const component )
{
  return field(index, component );
} 

template<>
GEOS_HOST_DEVICE
inline real64 bcFieldValue<arrayView1d<real64 const>>( arrayView1d<real64 const> const & field, 
                                                       localIndex const index, 
                                                       int const )
{
  return field[index];
}


template< typename T >
inline real64 fieldLinearIndex( T const & field, 
                                localIndex const index0,
                                localIndex const index1 )
{
  return field.linearIndex(index0, index1 );
} 

template<>
inline real64 fieldLinearIndex<arrayView1d<real64 const>>( arrayView1d<real64 const> const & field, 
                                                           localIndex const index0,
                                                           localIndex const index1 )
{
  return index0;
}


/**
 * @brief Wrapper for a linear operator applying the boundary conditions in a matrix-free fashion.
 * 
 * @tparam Vector The type of vectors used.
 * @tparam PrimaryFieldType The type of the primary field.
 */
template < typename Vector, typename PrimaryFieldType >
class LinearOperatorWithBC : public LinearOperator< Vector >
{
public:
  enum class DiagPolicy : integer
  {
    KeepDiagonal, ///< Use diagonal values
    DiagonalOne,  ///< Use one on the diagonal
    DiagonalZero, ///< Use zero on the diagonal
  };

  // struct fieldFunctor
  // {
  //   GEOS_HOST_DEVICE real64 operator()( localIndex const a )
  //   {
  //     return 
  //   }

  //   int const component;
  //   int const numComponents;
  // };

  /**
   * @brief Construct a new linear operator with boundary conditions wrapping the given linear operator.
   * 
   * @param solver 
   * @param unconstrained_op The linear operator without boundary conditions.
   * @param domain 
   * @param dofManager 
   * @param fieldName 
   * @param time 
   * @param diagPolicy The digonal policy to apply the boundary conditions.
   */
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

    // Compute the "diagonal"
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
        GEOS_ERROR("Not yet implemented.");
        // TODO: compute diagonal
        break;
    }

    // compute m_constrainedDofIndices and m_rhsContributions
    FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();
    globalIndex totalSize = 0;
    solver.forDiscretizationOnMeshTargets( m_domain.getMeshBodies(), [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
    {
      fsManager.apply<NodeManager>( m_time,
                       mesh,
                       m_fieldName,
                       [&]( FieldSpecificationBase const & bc,
                            string const &,
                            SortedArrayView< localIndex const > const & targetSet,
                            dataRepository::Group & targetGroup,
                            string const & GEOS_UNUSED_PARAM( fieldName ) )
      {
        totalSize += targetSet.size();
      } );
    } );

    
    // NOTE: we're not checking for duplicates.
    m_constrainedIndices.reserve( totalSize );
    m_constrainedDofIndices.reserve( totalSize );
    m_rhsContributions.reserve( totalSize );

    setupBoundaryConditions( solver, fsManager );
    // TODO: sort m_constrainedDofIndices and m_rhsContributions together

    srcWithBC.create( dofManager.numLocalDofs(), this->comm() );
    tmpRhs.create( dofManager.numLocalDofs(), this->comm() );
  }

  /**
   * @brief Compute the boundary condition.
   * 
   * @param solver 
   * @param fsManager 
   */
  void setupBoundaryConditions( SolverBase const & solver,
                                FieldSpecificationManager const & fsManager )
  {
    solver.forDiscretizationOnMeshTargets( m_domain.getMeshBodies(), [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
    {
      auto const & nodeManager = mesh.getNodeManager();
      auto const & field = nodeManager.getReference< PrimaryFieldType >( m_fieldName ).toViewConst();

      fsManager.apply<NodeManager>( m_time,
                       mesh,
                       m_fieldName,
                       [&]( FieldSpecificationBase const & bc,
                            string const &,
                            SortedArrayView< localIndex const > const & targetSet,
                            dataRepository::Group & targetGroup,
                            string const & fieldName )
      {
        array1d< localIndex > indexArray( targetSet.size() );
        int count = 0;
        for( localIndex const index : targetSet )
        {
          indexArray[ count++ ] = fieldLinearIndex( field, index, bc.getComponent() ) ;
        }


        array1d< globalIndex > dofArray( targetSet.size() );
        dofArray.setName("dofArray");
        arrayView1d< globalIndex > const & dof = dofArray.toView();

        array1d< real64 > rhsContributionArray( targetSet.size() );
        rhsContributionArray.setName("rhsContributionArray");
        arrayView1d< real64 > const & rhsContribution = rhsContributionArray.toView(); 
        arrayView1d< globalIndex const > const & dofMap = targetGroup.getReference< array1d< globalIndex > >( m_dofManager.getKey( fieldName ) );
        int const component = bc.getComponent();
        int const rankOffset = m_dofManager.rankOffset();
        // we need to write a new variant of this function that calculates the rhs contribution to the original array2d instead of the dof vector.
        bc.computeRhsContribution< FieldSpecificationEqual, 
                                   parallelDevicePolicy<> >( targetSet,
                                                             m_time,
                                                             1.0, // TODO: double check
                                                             targetGroup,
                                                             dofMap,
                                                             rankOffset,
                                                             m_diagonal.toViewConst(),
                                                             dof,
                                                             rhsContribution,
                                                             [=] GEOS_HOST_DEVICE (localIndex const a)->real64
                                                             {
                                                               return bcFieldValue( field, a, component );
                                                             } );
                                                             //field );

        dof.move( LvArray::MemorySpace::host, false );
        rhsContribution.move( LvArray::MemorySpace::host, false );

        m_constrainedDofIndices.insert( m_constrainedDofIndices.size(), dofArray.begin(), dofArray.end() );
        m_rhsContributions.insert( m_rhsContributions.size(), rhsContribution.begin(), rhsContribution.end() );
        m_constrainedIndices.insert( m_constrainedIndices.size(), indexArray.begin(), indexArray.end() );
        
      } );
    } );
  }

  /**
   * @brief Compute the contributions of the boundary conditions on the right hand side and the
   * initial guess.
   * 
   * @param rhs The right hand side vector.
   * @param solution The initial guess vector.
   */
  void computeConstrainedRHS( ParallelVector & rhs, ParallelVector & solution ) const
  {
    GEOS_MARK_FUNCTION;

    using POLICY = parallelDevicePolicy<>;

    // Construct [x_BC,0]
    srcWithBC.zero();
    arrayView1d< real64 > const localBC = srcWithBC.open();
    arrayView1d< real64 > const initSolution = solution.open();
    arrayView1d< localIndex const > const localBCIndices = m_constrainedIndices.toViewConst();
    arrayView1d< localIndex const > const localBCDofs = m_constrainedDofIndices.toViewConst();
    arrayView1d< real64 const > const localDiag = m_diagonal.toViewConst();
    arrayView1d< real64 const > const localRhsContributions = m_rhsContributions.toViewConst();
    forAll< POLICY >( m_constrainedDofIndices.size(),
                      [ initSolution, localBC, localBCIndices, localBCDofs, localDiag, localRhsContributions ] GEOS_HOST_DEVICE
                        ( localIndex const i )
    {
      localIndex const idx = localBCIndices[ i ];
      localBC[ idx ] = localRhsContributions [ i ] / localDiag[ localBCDofs[i] ];
      initSolution[idx] = localBC[ idx ];
    } );
    srcWithBC.close();
    solution.close();


    // Bottom contribution to rhs
    tmpRhs = rhs;
    
    // std::cout<<"LinearOperatorWithBC::computeConstrainedRHS - srcWithBC"<<std::endl;
    // std::cout<<srcWithBC<<std::endl;
    m_unconstrained_op.apply( srcWithBC, rhs );

    // std::cout<<"LinearOperatorWithBC::computeConstrainedRHS - rhs"<<std::endl;
    // std::cout<<rhs<<std::endl;


    rhs.axpby( 1.0, tmpRhs, -1.0 );

    // D_GG x_BC
    arrayView1d< real64 > const localRhs = rhs.open();
    forAll< POLICY >( m_constrainedDofIndices.size(),
                      [ localRhs, localBCIndices, localRhsContributions ] GEOS_HOST_DEVICE
                        ( localIndex const i )
    {
      localIndex const idx = localBCIndices[ i ]; 
      localRhs[ idx ] = localRhsContributions [ i ];
    } );
    rhs.close();
  }

  /**
   * @brief Apply the linear operator with the boundary conditions.
   * 
   * @param src The input vector.
   * @param dst The output vector.
   */
  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const
  {
    GEOS_MARK_FUNCTION;

    using POLICY = parallelDeviceAsyncPolicy<>;

    arrayView1d< real64 > const localSrcWithBC = srcWithBC.open();
    arrayView1d< real64 const > const localSrc = src.values();

    forAll< POLICY >( localSrc.size(), [localSrcWithBC, localSrc] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localSrcWithBC[ i ] = localSrc[ i ];
    } );

    arrayView1d< localIndex const > const localBCIndices = m_constrainedIndices.toViewConst();
//    std::cout << "constrained ind: " << localBCIndices << std::endl;

    forAll< POLICY >( m_constrainedDofIndices.size(), [localSrcWithBC, localBCIndices] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localSrcWithBC[ localBCIndices[ i ] ] = 0.0;
    } );
    
    srcWithBC.close();

//    std::cout << "srcWithBC: " << srcWithBC << std::endl;

    m_unconstrained_op.apply( srcWithBC, dst );

//    std::cout << "dst: " << dst << std::endl;

//    std::cout << "localSrc: " << localSrc << std::endl;
    arrayView1d< real64 > const localDst = dst.open();
    switch (m_diagPolicy)
    {
      case DiagPolicy::DiagonalOne:
        {
          forAll< POLICY >( m_constrainedDofIndices.size(), [ localSrc, localDst, localBCIndices ] GEOS_HOST_DEVICE ( localIndex const i )
          {
            localIndex const idx = localBCIndices[ i ];
            localDst[ idx ] = localSrc [ idx ];
          } );
//          std::cout << "localDst: " << localDst << std::endl;
        }
        break;
      case DiagPolicy::DiagonalZero:
        {
          forAll< POLICY >( m_constrainedDofIndices.size(), [ localDst, localBCIndices ] GEOS_HOST_DEVICE ( localIndex const i )
          {
            localIndex const idx = localBCIndices[ i ];
            localDst[ idx ] = 0.0;
          } );
        }
        break;
      case DiagPolicy::KeepDiagonal:
        GEOS_ERROR("Not yet implemented.");
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
  array1d<localIndex> m_constrainedIndices;
  array1d<real64> m_rhsContributions;
  array1d< localIndex > m_constrainedDofIndices;

  mutable ParallelVector srcWithBC;
  mutable ParallelVector tmpRhs;
  array1d< real64 > m_diagonal;
};

} /* namespace geos */

#endif /* GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATORWITHBC_HPP_ */
