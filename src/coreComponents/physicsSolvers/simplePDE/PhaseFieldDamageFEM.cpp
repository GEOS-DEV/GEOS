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
 * @file PhaseFieldDamageFEM.cpp
 *
 */

#include "PhaseFieldDamageFEM.hpp"

#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"

#include "managers/DomainPartition.hpp"

// this should be part of the input file

double myFunc2( double, double, double )
{
  return 0;
  // return pow(x, 2) + pow(y, 2) + pow(z, 2) + 6;
//  return x * (1 - x) * y * (1 - y) * z * (1 - z) -
//         2 * (x - 1) * x * (y - 1) * y - 2 * (x - 1) * x * (z - 1) * z -
//         2 * (y - 1) * y * (z - 1) * z;
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

PhaseFieldDamageFEM::PhaseFieldDamageFEM( const std::string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_solidModelName()
{

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.Key() )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "option for default time integration method" );

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.Key(), &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "name of field variable" );

  registerWrapper( viewKeyStruct::localDissipationOption, &m_localDissipationOption )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription("Type of local dissipation function. Can be Linear or Quadratic" );

  registerWrapper(viewKeyStruct::lengthScale, &m_lengthScale)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("lenght scale l in the phase-field equation");

  registerWrapper(viewKeyStruct::criticalFractureEnergy, &m_criticalFractureEnergy)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("critical fracture energy");

  registerWrapper< string >( viewKeyStruct::solidModelNameString, &m_solidModelName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "name of solid constitutive model" );
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel *meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    NodeManager * const nodes = meshLevel->getNodeManager();

    nodes->registerWrapper< real64_array >( m_fieldName )
      ->setApplyDefaultValue( 0.0 )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setDescription( "Primary field variable" );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion >( [ &]( CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper( viewKeyStruct::coeffName, &m_coeff )->
        setApplyDefaultValue( 0.0 )->
        setPlotLevel( PlotLevel::LEVEL_0 )->
        setDescription( "field variable representing the diffusion coefficient" );
    } );

  }
}

void PhaseFieldDamageFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  string tiOption = this->getReference< string >(
    PhaseFieldDamageFEMViewKeys.timeIntegrationOption );

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

  if( m_localDissipationOption != "Linear" && m_localDissipationOption != "Quadratic" )
  {
    GEOSX_ERROR( "invalid local dissipation option - must be Linear or Quadratic" );
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

real64 PhaseFieldDamageFEM::SolverStep( real64 const & time_n, real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition *domain )
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
      this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager,
                                   m_matrix,
                                   m_rhs, m_solution );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::ExplicitStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::ImplicitStepSetup(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition * const domain,
  DofManager & dofManager,
  ParallelMatrix & matrix,
  ParallelVector & rhs, ParallelVector & solution )
{
  // Computation of the sparsity pattern
  SetupSystem( domain, dofManager, matrix, rhs, solution );
}

void PhaseFieldDamageFEM::ImplicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::SetupDofs(
  DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName, DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );

}

void PhaseFieldDamageFEM::AssembleSystem( real64 const time_n,
                                          real64 const GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );

  arrayView1d< globalIndex const > const & dofIndex = nodeManager->getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );

  //arrayView1d<R1Tensor> &X = nodeManager->referencePosition();

  // Initialize all entries to zero
  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  // begin region loop
  for( localIndex er = 0; er < elemManager->numRegions(); ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion( er );

    elementRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( esr ),
                                                                           CellElementSubRegion & elementSubRegion )
    {

      constitutive::ConstitutiveBase * const
      solidModel = elementSubRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( m_solidModelName );

      constitutive::ConstitutivePassThru< constitutive::DamageBase >::Execute( solidModel,
                                                                               [&]( auto * const damageModel )
      {
        using CONSTITUTIVE_TYPE = TYPEOFPTR( damageModel );
        typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel->createKernelUpdates();

        arrayView4d< real64 const > const &
        dNdX = elementSubRegion.dNdX();

        arrayView2d< real64 const > const &
        detJ = elementSubRegion.detJ();

        localIndex const numNodesPerElement =  elementSubRegion.numNodesPerElement();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemNodes = elementSubRegion.nodeList();

        // arrayView1d<real64 const> const &
        // coeff = elementSubRegion.getReference<array1d<real64> >(viewKeyStruct::coeffName);

        globalIndex_array elemDofIndex( numNodesPerElement );
        real64_array element_rhs( numNodesPerElement );
        real64_array2d element_matrix( numNodesPerElement, numNodesPerElement );

        arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();
        std::unique_ptr< FiniteElementBase > finiteElement = feDiscretization->getFiniteElement( elementSubRegion.GetElementTypeString() );
        localIndex const n_q_points = finiteElement->n_quadrature_points();

        real64 ell = m_lengthScale;                       //phase-field length scale
        real64 Gc = m_criticalFractureEnergy;             //energy release rate
        double threshold = 3 * Gc / (16 * ell);           //elastic energy threshold - use when Local Dissipation is linear

        arrayView1d< real64 > const & nodalDamage = nodeManager->getReference< array1d< real64 > >( m_fieldName );
        //real64 diffusion = 1.0;
        // begin element loop, skipping ghost elements
        for( localIndex k = 0; k < elementSubRegion.size(); ++k )
        {
          if( elemGhostRank[k] < 0 )
          {
            element_rhs = 0.0;
            element_matrix = 0.0;
            for( localIndex q = 0; q < n_q_points; ++q )
            {
              real64 const strainEnergyDensity = constitutiveUpdate.calculateStrainEnergyDensity( k,q );
              real64 D = 0;                                                                   //max between threshold and
                                                                                              // Elastic energy
              if( m_localDissipationOption == "Linear" )
              {
                D = std::max( threshold, strainEnergyDensity );
                //D = max(strainEnergy(k,q), strainEnergy(k,q));//debbuging line - remove after testing
              }
              //Interpolate d and grad_d

              real64 qp_damage = 0.0;
              R1Tensor qp_grad_damage;
              R1Tensor temp;
              for( localIndex a = 0; a < numNodesPerElement; ++a )
              {
                qp_damage += finiteElement->value( a, q ) * nodalDamage[elemNodes( k, a )];
                temp = dNdX[k][q][a];
                temp *= nodalDamage[elemNodes( k, a )];
                qp_grad_damage += temp;

              }
              //std::cout << "Damage: " << qp_damage <<std::endl;
              //std::cout << "GradDamage: " << qp_grad_damage <<std::endl;
              for( localIndex a = 0; a < numNodesPerElement; ++a )
              {
                elemDofIndex[a] = dofIndex[elemNodes( k, a )];
                //real64 diffusion = 1.0;
                real64 Na = finiteElement->value( a, q );
                //element_rhs(a) += detJ[k][q] * Na * myFunc(Xq, Yq, Zq); //older reaction diffusion solver
                if( m_localDissipationOption == "Linear" )
                {
                  element_rhs( a ) += detJ[k][q] * (Na * (ell * D - 3 * Gc / 16 )/ Gc -
                                                    0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( qp_grad_damage, dNdX[k][q][a] ) -
                                                    (ell * D/Gc) * Na * qp_damage);
                }
                else
                {
                  element_rhs( a ) += detJ[k][q] * (Na * (2 * ell) * strainEnergyDensity / Gc -
                                                    (pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( qp_grad_damage, dNdX[k][q][a] ) +
                                                     Na * qp_damage * (1 + 2 * ell*strainEnergyDensity/Gc)) );
                }
                for( localIndex b = 0; b < numNodesPerElement; ++b )
                {
                  real64 Nb = finiteElement->value( b, q );
                  if( m_localDissipationOption == "Linear" )
                  {
                    element_matrix( a, b ) -= detJ[k][q] *
                                              (0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( dNdX[k][q][a], dNdX[k][q][b] ) +
                                               (ell * D/Gc) * Na * Nb);
                  }
                  else
                  {
                    element_matrix( a, b ) -= detJ[k][q] *
                                              ( pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( dNdX[k][q][a], dNdX[k][q][b] ) +
                                                  Na * Nb * (1 + 2 * ell*strainEnergyDensity/Gc )
                                              );
                  }
                }
              }
            }
            matrix.add( elemDofIndex, elemDofIndex, element_matrix );
            rhs.add( elemDofIndex, element_rhs );
          }
        }
      } );
    } );
  }
  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::AssembleSystem" );
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

    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void PhaseFieldDamageFEM::ApplySystemSolution( DofManager const & dofManager,
                                               ParallelVector const & solution,
                                               real64 const scalingFactor,
                                               DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( solution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );

  // Syncronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getNeighbors() );
}

void PhaseFieldDamageFEM::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition * const domain,
  DofManager const & dofManager,
  ParallelMatrix & matrix, ParallelVector & rhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, *domain, m_matrix, m_rhs );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
    integer newtonIter = solverParams.m_numNewtonIterations;

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" +
                          std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" +
                          std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs );

    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

real64
PhaseFieldDamageFEM::CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                            DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                            ParallelVector const & rhs )
{

  real64 const norm = rhs.norm2();
  return norm;

}

void PhaseFieldDamageFEM::SolveSystem( DofManager const & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::SolveSystem" );
    GEOSX_LOG_RANK_0( "\nSolution\n" );
    std::cout << solution;
  }
}

void PhaseFieldDamageFEM::ApplyDirichletBC_implicit( real64 const time,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     ParallelMatrix & matrix,
                                                     ParallelVector & rhs )

{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  matrix.open();
  rhs.open();
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

  fsManager.ApplyFieldValue< serialPolicy >( time, &domain, "ElementRegions", viewKeyStruct::coeffName );
  rhs.close();
  matrix.close();
}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldDamageFEM, std::string const &,
                        Group * const )
} // namespace geosx
