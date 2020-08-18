/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior
 * University Copyright (c) 2018-2019 Total, S.A Copyright (c) 2019-     GEOSX
 * Contributors All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS
 * files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ReactionDiffusionFEM.cpp
 *
 */

#include "ReactionDiffusionFEM.hpp"

#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"

#include "managers/DomainPartition.hpp"

// this should be part of the input file

double myFunc( double x, double y, double z )
{
  //return 0;
  // return pow(x, 2) + pow(y, 2) + pow(z, 2) + 6;
  int N = 8;
  return x * y * z * ( N - x ) * ( N - y ) * ( N - z ) - 2 * y * z * ( N - y ) * ( N - z ) - 2 * x * y * ( N - x ) * ( N - y ) - 2 * x * z * ( N - x ) *
         ( N - z );
  // return (1 - x*x) * (1 - y*y) * (1 - z*z) -
  //        2 * (1 - x*x) * (1 - y*y) - 2 * (1 - x*x) * (1 - z*z) -
  //        2 * (1 - y*y) * (1 - z*z);
}

///////////////////////////////////

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
} // namespace dataRepository

using namespace dataRepository;
using namespace constitutive;

ReactionDiffusionFEM::ReactionDiffusionFEM( const std::string & name,
                                            Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" )
{

  registerWrapper< string >( reactionDiffusionFEMViewKeys.timeIntegrationOption.Key() )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "option for default time integration method" );

  registerWrapper< string >( reactionDiffusionFEMViewKeys.fieldVarName.Key(), &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "name of field variable" );

//  registerWrapper(viewKeyStruct::coeffFieldName,
//                          &m_coeffFieldName, false)->
//    setInputFlag(InputFlags::REQUIRED)->
//    setDescription("name of field variable representing the diffusion coefficient");

}

ReactionDiffusionFEM::~ReactionDiffusionFEM()
{
  // TODO Auto-generated destructor stub
}

void ReactionDiffusionFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    //MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    NodeManager * const nodes = mesh.second->group_cast< MeshBody * >()
                                  ->getMeshLevel( 0 )
                                  ->getNodeManager();

    nodes->registerWrapper< real64_array >( m_fieldName )
      ->setApplyDefaultValue( 0.0 )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setDescription( "Primary field variable" );
    //use this to get coeff from input file
    //ElementRegionManager * const elemManager = meshLevel->getElemManager();

    // elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    // {
    //   subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::coeffName)->
    //     setApplyDefaultValue(0.0)->
    //     setPlotLevel(PlotLevel::LEVEL_0)->
    //     setDescription("field variable representing the diffusion coefficient");
    // });

  }
}

void ReactionDiffusionFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  string tiOption = this->getReference< string >(
    reactionDiffusionFEMViewKeys.timeIntegrationOption );

  if( tiOption == "SteadyState" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
  }
  else if( tiOption == "ImplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
  }
  else if( tiOption == "ExplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
  }
  else
  {
    GEOSX_ERROR( "invalid time integration option" );
  }

  // Set basic parameters for solver
  // m_linearSolverParameters.logLevel = 0;
  // m_linearSolverParameters.solverType = "gmres";
  // m_linearSolverParameters.krylov.tolerance = 1e-8;
  // m_linearSolverParameters.krylov.maxIterations = 250;
  // m_linearSolverParameters.krylov.maxRestart = 250;
  // m_linearSolverParameters.preconditionerType = "amg";
  // m_linearSolverParameters.amg.smootherType = "gaussSeidel";
  // m_linearSolverParameters.amg.coarseType = "direct";
}

real64 ReactionDiffusionFEM::SolverStep( real64 const & time_n, real64 const & dt,
                                         const int cycleNumber,
                                         DomainPartition & domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption ==
           timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    dtReturn =
      this->LinearImplicitStep( time_n,
                                dt,
                                cycleNumber,
                                domain );
  }
  return dtReturn;
}

real64 ReactionDiffusionFEM::ExplicitStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}

void ReactionDiffusionFEM::ImplicitStepSetup(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & domain)
{
  // Computation of the sparsity pattern
  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void ReactionDiffusionFEM::ImplicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ReactionDiffusionFEM::SetupDofs(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName, DofManager::Location::Node );
}

void ReactionDiffusionFEM::AssembleSystem( real64 const time_n,
                                           real64 const GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & domain,
                                           DofManager const & dofManager,
                                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs)
{
  MeshLevel * const mesh = domain.getMeshBody(0)->getMeshLevel(0);

  NodeManager * const nodeManager = mesh->getNodeManager();

  ElementRegionManager * const elemManager = mesh->getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );

  arrayView1d< globalIndex const > const & dofIndex =
    nodeManager->getReference< array1d< globalIndex > >(
      dofManager.getKey( m_fieldName ) );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();


  // begin region loop
  for( localIndex er = 0; er < elemManager->numRegions(); ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion( er );

    elementRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( esr ),
                                                                           CellElementSubRegion const & elementSubRegion )
    {
      arrayView4d< real64 const > const & dNdX = elementSubRegion.dNdX();

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

      localIndex const numNodesPerElement =  elementSubRegion.numNodesPerElement();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const &
      elemNodes = elementSubRegion.nodeList();

      //use this to get coeff from input file
      //arrayView1d<real64 const> const &
      //coeff = elementSubRegion.getReference<array1d<real64> >(viewKeyStruct::coeffName);

      globalIndex_array elemDofIndex( numNodesPerElement );
      real64_array element_rhs( numNodesPerElement );
      real64_array2d element_matrix( numNodesPerElement, numNodesPerElement );

      arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();

      std::unique_ptr< FiniteElementBase > finiteElement = feDiscretization->getFiniteElement( elementSubRegion.GetElementTypeString() );
      localIndex const n_q_points = finiteElement->n_quadrature_points();

      real64 reaction = 1.0;
      real64 diffusion = 1.0;
      // begin element loop, skipping ghost elements
      for( localIndex k = 0; k < elementSubRegion.size(); ++k )
      {
        if( elemGhostRank[k] < 0 )
        {
          element_rhs = 0.0;
          element_matrix = 0.0;
          for( localIndex q = 0; q < n_q_points; ++q )
          {
            real64 Xq = 0;
            real64 Yq = 0;
            real64 Zq = 0;
            for( localIndex a = 0; a < numNodesPerElement; ++a )
            {
              Xq = Xq + finiteElement->value( a, q ) * X[elemNodes( k, a )][0];

              Yq = Yq + finiteElement->value( a, q ) * X[elemNodes( k, a )][1];

              Zq = Zq + finiteElement->value( a, q ) * X[elemNodes( k, a )][2];
            }
            for( localIndex a = 0; a < numNodesPerElement; ++a )
            {
              elemDofIndex[a] = dofIndex[elemNodes( k, a )];

              real64 Na = finiteElement->value( a, q );
              //may need a minus sign here
              element_rhs( a ) += -detJ[k][q] * Na * myFunc( Xq, Yq, Zq );                                                              //older
                                                                                                                                        // reaction
                                                                                                                                        // diffusion
                                                                                                                                        // solver
              for( localIndex b = 0; b < numNodesPerElement; ++b )
              {
                real64 Nb = finiteElement->value( b, q );
                element_matrix( a, b ) +=
                  detJ[k][q] *
                  (diffusion * - LvArray::tensorOps::AiBi<3>( dNdX[k][q][a], dNdX[k][q][b] ) + reaction * Na * Nb);

              }
            }
          }
          matrix.add( elemDofIndex, elemDofIndex, element_matrix );
          rhs.add( elemDofIndex, element_rhs );
        }
      }
    } );
  }
  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After ReactionDiffusionFEM::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
    integer newtonIter = solverParams.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" +
                          std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" +
                          std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs );

    GEOSX_LOG_RANK_0( "After ReactionDiffusionFEM::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void ReactionDiffusionFEM::ApplySystemSolution( DofManager const & dofManager,
                                                arrayView1d< real64 const > const & localSolution,
                                                real64 const scalingFactor,
                                                DomainPartition & domain )
{
  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.copyVectorToField( localSolution,
                                m_fieldName,
                                m_fieldName,
                                scalingFactor );

  // Syncronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  CommunicationTools::SynchronizeFields(
    fieldNames,
    mesh,
    domain.getNeighbors() );
}

void ReactionDiffusionFEM::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex > const & localMatrix,
                    arrayView1d< real64 > const & localRhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, domain, localMatrix, localRhs );

//  if( getLogLevel() == 2 )
//  {
//    GEOSX_LOG_RANK_0( "After ReactionDiffusionFEM::ApplyBoundaryConditions" );
//    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
//    std::cout << matrix;
//    GEOSX_LOG_RANK_0( "\nResidual:\n" );
//    std::cout << rhs;
//  }
//
//  if( getLogLevel() >= 3 )
//  {
//    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
//    integer newtonIter = solverParams.m_numNewtonIterations;
//
//    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    matrix.write( filename_mat );
//
//    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    rhs.write( filename_rhs );
//
//    GEOSX_LOG_RANK_0( "After ReactionDiffusionFEM::ApplyBoundaryConditions" );
//    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
//    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
//  }
}

void ReactionDiffusionFEM::SolveSystem( DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs,
                                        ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After ReactionDiffusionFEM::SolveSystem" );
    GEOSX_LOG_RANK_0( "\nSolution\n" );
    std::cout << solution;
  }
}

void ReactionDiffusionFEM::ApplyDirichletBC_implicit(
  real64 const time,
  DofManager const & dofManager, DomainPartition & domain,
  ParallelMatrix & matrix,
  ParallelVector & rhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const * const bc, string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const GEOSX_UNUSED_PARAM( fieldName ) ) -> void
  {
    bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                        LAInterface >( targetSet,
                                                       time,
                                                       targetGroup,
                                                       m_fieldName,
                                                       dofManager.getKey( m_fieldName ),
                                                       1,
                                                       matrix,
                                                       rhs );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, ReactionDiffusionFEM, std::string const &,
                        Group * const )
} // namespace geosx
