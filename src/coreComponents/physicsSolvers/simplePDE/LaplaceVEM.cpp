/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LaplaceVEM.cpp
 *
 */

#include "LaplaceVEM.hpp"

#include <vector>
#include <math.h>

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "dataRepository/Group.hpp"
#include "common/TimingMacros.hpp"

#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"
#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

namespace geosx
{

  namespace dataRepository
  {
    namespace keys
    {}
  }

  using namespace dataRepository;
  using namespace constitutive;


  //START_SPHINX_INCLUDE_01
  LaplaceVEM::LaplaceVEM( const std::string& name,
  Group * const parent ):
    SolverBase( name, parent ),
    m_fieldName("primaryField")
  {
    registerWrapper<string>(laplaceVEMViewKeys.timeIntegrationOption.Key())->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("option for default time integration method");

    registerWrapper<string>(laplaceVEMViewKeys.fieldVarName.Key(), &m_fieldName, false)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("name of field variable");
  }
  //END_SPHINX_INCLUDE_01

  LaplaceVEM::~LaplaceVEM()
  {
    // TODO Auto-generated destructor stub
  }


  //START_SPHINX_INCLUDE_02
  void LaplaceVEM::RegisterDataOnMesh( Group * const MeshBodies )
  {
    for( auto & mesh : MeshBodies->GetSubGroups() )
    {
      NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();

      nodes->registerWrapper<real64_array >( m_fieldName )->
        setApplyDefaultValue(0.0)->
        setPlotLevel(PlotLevel::LEVEL_0)->
        setDescription("Primary field variable");
    }
  }
  //END_SPHINX_INCLUDE_02

  //START_SPHINX_INCLUDE_03
  void LaplaceVEM::PostProcessInput()
  {
    SolverBase::PostProcessInput();

    string tiOption = this->getReference<string>(laplaceVEMViewKeys.timeIntegrationOption);

    if( tiOption == "SteadyState" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
    }
    else if( tiOption == "ImplicitTransient" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
    }
    else if ( tiOption == "ExplicitTransient" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
    }
    else
    {
      GEOSX_ERROR("invalid time integration option");
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
  //END_SPHINX_INCLUDE_03

  real64 LaplaceVEM::SolverStep( real64 const& time_n,
  real64 const& dt,
  const int cycleNumber,
  DomainPartition * domain )
  {
    real64 dtReturn = dt;
    if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
    {
      dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );
    }
    else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitTransient ||
    m_timeIntegrationOption == timeIntegrationOption::SteadyState )
    {
      dtReturn = this->LinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );
    }
    return dtReturn;
  }

  real64 LaplaceVEM::ExplicitStep( real64 const& GEOSX_UNUSED_ARG( time_n ),
  real64 const& dt,
  const int GEOSX_UNUSED_ARG( cycleNumber ),
  DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
  {
    return dt;
  }

  void LaplaceVEM::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
  real64 const & GEOSX_UNUSED_ARG( dt ),
  DomainPartition * const domain,
  DofManager & dofManager,
  ParallelMatrix & matrix,
  ParallelVector & rhs,
  ParallelVector & solution )
  {
    // Computation of the sparsity pattern
    SetupSystem( domain, dofManager, matrix, rhs, solution );
  }

  void LaplaceVEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG( time_n ),
  real64 const & GEOSX_UNUSED_ARG( dt ),
  DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
  {
  }

  void LaplaceVEM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
  DofManager & dofManager ) const
  {
    dofManager.addField( m_fieldName,
    DofManager::Location::Node// ,
    // DofManager::Connectivity::Elem
                         );
  }

  //START_SPHINX_INCLUDE_04
  void LaplaceVEM::AssembleSystem( real64 const time_n,
  real64 const GEOSX_UNUSED_ARG( dt ),
  DomainPartition * const domain,
  DofManager const & dofManager,
  ParallelMatrix & matrix,
  ParallelVector & rhs )
  {
    MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
    NodeManager const * const nodeManager = mesh->getNodeManager();
    ElementRegionManager * const elemManager = mesh->getElemManager();
    NumericalMethodsManager const *
      numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
    FiniteElementDiscretizationManager const *
      feDiscretizationManager = numericalMethodManager->
      GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

    array1d<globalIndex> const & dofIndex =
      nodeManager->getReference< array1d<globalIndex> >( dofManager.getKey( m_fieldName ) );

    // Initialize all entries to zero
    matrix.zero();
    rhs.zero();

    matrix.open();
    rhs.open();

    // get node properties
    arrayView1d<R1Tensor const> const & nodeCoords = nodeManager->referencePosition();

    // get face properties
    FaceManager const * const faceManager = mesh->getFaceManager();
    arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
    GEOSX_UNUSED_VAR(faceArea);
    arrayView1d<R1Tensor const> const & faceCenter = faceManager->faceCenter();
    GEOSX_UNUSED_VAR(faceCenter);
    arrayView1d<R1Tensor const> const & faceNormals = faceManager->faceNormal();

    // get edge properties
    EdgeManager const * const edgeManager = mesh->getEdgeManager();
    R1Tensor edgeCenter;
    edgeManager->calculateCenter(0, nodeCoords, edgeCenter);
    GEOSX_UNUSED_VAR(edgeCenter);

    // begin region loop
    for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
    {
      ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

      FiniteElementDiscretization const *
        feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

      elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const GEOSX_UNUSED_ARG( esr ),
      CellElementSubRegion const * const elementSubRegion )
      {
        array3d<R1Tensor> const &
          dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

        arrayView2d<real64> const &
          detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

        localIndex const numNodesPerElement = elementSubRegion->numNodesPerElement();
        arrayView2d<localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM> const & elemNodes = elementSubRegion->nodeList();

        globalIndex_array elemDofIndex( numNodesPerElement );
        real64_array element_rhs( numNodesPerElement );
        real64_array2d element_matrix( numNodesPerElement, numNodesPerElement );

        integer_array const & elemGhostRank = elementSubRegion->m_ghostRank;
        localIndex const n_q_points = feDiscretization->m_finiteElement->n_quadrature_points();

        arrayView2d<localIndex> const & elemsToFaces = elementSubRegion->faceList(); GEOSX_UNUSED_VAR(elemsToFaces);
        ArrayOfArraysView< localIndex const > const & facesToEdges = faceManager->edgeList(); GEOSX_UNUSED_VAR(facesToEdges);
        ArrayOfArraysView< localIndex const > const & facesToNodes = faceManager->nodeList();
        arrayView2d<localIndex> const & edgesToNodes = edgeManager->nodeList(); GEOSX_UNUSED_VAR(edgesToNodes);


        localIndex const & numFacesPerElem = elemsToFaces.size(1);

        // begin element loop, skipping ghost elements
        for( localIndex k=0 ; k<elementSubRegion->size() ; ++k )
        {
          // compute basis functions integral on each face
          real64_array2d basisFunctionsIntegrals( numFacesPerElem, numNodesPerElement ); GEOSX_UNUSED_VAR(basisFunctionsIntegrals);
          for( localIndex kf = 0; kf < numFacesPerElem; ++kf )
          {
            localIndex const & faceIndex = elemsToFaces(k,kf);
            localIndex const & numFaceNodes = facesToNodes.sizeOfArray(faceIndex);
            // compute face rotation matrix
            R1TensorT<3> const & V0 = nodeCoords(facesToNodes(faceIndex, 0));
            R1TensorT<3> V1mV0 = nodeCoords(facesToNodes(faceIndex, 1));
            V1mV0 -= V0;
            realT normV1mV0 = V1mV0.L2_Norm();
            real64_array2d auxMatZ(numFaceNodes, 3);
            auxMatZ(0, 0) = V1mV0(0);
            auxMatZ(0, 1) = V1mV0(1);
            auxMatZ(0, 2) = V1mV0(2);
            real64_array2d auxMatW(3, numFaceNodes);
            auxMatW(0, 0) = normV1mV0;
            auxMatW(1, 0) = 0;
            auxMatW(2, 0) = 0;
            for( localIndex col = 1; col < numFaceNodes-1; ++col )
            {
              R1TensorT<3> VimV0 = nodeCoords(facesToNodes(faceIndex, col+1));
              VimV0 -= V0;
              realT normVimV0 = VimV0.L2_Norm(); GEOSX_UNUSED_VAR(normVimV0);
              auxMatZ(col, 0) = VimV0(0);
              auxMatZ(col, 1) = VimV0(1);
              auxMatZ(col, 2) = VimV0(2);
              realT inverseProdNorms = 1.0/(normV1mV0*normVimV0);
              realT cosAngleBetweenVectors = Dot(VimV0,V1mV0)*inverseProdNorms; GEOSX_UNUSED_VAR(cosAngleBetweenVectors);
              realT sinAngleBetweenVectors = Cross(VimV0,V1mV0).L2_Norm()*inverseProdNorms;
              auxMatW(0, col) = normVimV0*cosAngleBetweenVectors;
              auxMatW(1, col) = normVimV0*sinAngleBetweenVectors;
              auxMatW(2, col) = 0;
            }
            R1TensorT<3> const & faceNormal = faceNormals(faceIndex); GEOSX_UNUSED_VAR(faceNormal);
            auxMatZ(numFaceNodes-1, 0) = faceNormal(0);
            auxMatZ(numFaceNodes-1, 1) = faceNormal(1);
            auxMatZ(numFaceNodes-1, 2) = faceNormal(2);
            auxMatW(0, numFaceNodes-1) = 0;
            auxMatW(1, numFaceNodes-1) = 0;
            auxMatW(2, numFaceNodes-1) = 1;
            array2d<real64> auxMatH(3,3);
            for(localIndex i = 0; i < 9; ++i)
              auxMatH.data()[i] = 0;
            for(localIndex i = 0; i < numFaceNodes; ++i)
            {
              auxMatH(0,0) += auxMatW(0,i)*auxMatZ(i,0);
              auxMatH(0,1) += auxMatW(0,i)*auxMatZ(i,1);
              auxMatH(0,2) += auxMatW(0,i)*auxMatZ(i,2);
              auxMatH(1,0) += auxMatW(1,i)*auxMatZ(i,0);
              auxMatH(1,1) += auxMatW(1,i)*auxMatZ(i,1);
              auxMatH(1,2) += auxMatW(1,i)*auxMatZ(i,2);
              auxMatH(2,0) += auxMatW(2,i)*auxMatZ(i,0);
              auxMatH(2,1) += auxMatW(2,i)*auxMatZ(i,1);
              auxMatH(2,2) += auxMatW(2,i)*auxMatZ(i,2);
            }
            array2d<real64> svdMatU(3,3), svdMatVT(3,3);
            array1d<real64> svdValues(3);
            BlasLapackLA::matrixSVD(auxMatH, svdMatU, svdValues, svdMatVT);
            R2TensorT<3> faceRotationMatrix(0.0);
            // for(localIndex i = 0; i < 3; ++i)
            //   for(localIndex j = 0; j < 3; ++j)
          }


          GEOSX_UNUSED_VAR(dNdX);
          GEOSX_UNUSED_VAR(detJ);
          GEOSX_UNUSED_VAR(elemNodes);
          GEOSX_UNUSED_VAR(elemGhostRank);
          GEOSX_UNUSED_VAR(n_q_points);
          GEOSX_UNUSED_VAR(dofIndex);
          // if(elemGhostRank[k] < 0)
          // {
          //   element_rhs = 0.0;
          //   element_matrix = 0.0;
          //   for( localIndex q=0 ; q<n_q_points ; ++q)
          //   {
          //     for( localIndex a=0 ; a<numNodesPerElement ; ++a)
          //     {
          //       elemDofIndex[a] = dofIndex[ elemNodes( k, a ) ];

          //       real64 diffusion = 1.0;
          //       for( localIndex b=0 ; b<numNodesPerElement ; ++b)
          //       {
          //         element_matrix(a,b) += detJ[k][q] *
          //           diffusion *
          //           + Dot( dNdX[k][q][a], dNdX[k][q][b] );
          //       }

          //     }
          //   }
          //   matrix.add( elemDofIndex, elemDofIndex, element_matrix );
          //   rhs.add( elemDofIndex, element_rhs );
          // }
        }
      });
    }
    matrix.close();
    rhs.close();
    //END_SPHINX_INCLUDE_04

    // Debug for logLevel >= 2
    GEOSX_LOG_LEVEL_RANK_0( 2, "After LaplaceVEM::AssembleSystem" );
    GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
    GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );

    if( getLogLevel() >= 3 )
    {
      // SystemSolverParameters * const solverParams = getSystemSolverParameters();
      // integer newtonIter = solverParams->numNewtonIterations();
      integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

      string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
      matrix.write( filename_mat, true );

      string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
      rhs.write( filename_rhs, true );

      GEOSX_LOG_RANK_0( "After LaplaceVEM::AssembleSystem" );
      GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
      GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
    }
  }

  void LaplaceVEM::ApplySystemSolution( DofManager const & dofManager,
  ParallelVector const & solution,
  real64 const scalingFactor,
  DomainPartition * const domain )
  {
    dofManager.addVectorToField( solution, m_fieldName, m_fieldName, scalingFactor );

    // Synchronize ghost nodes
    std::map<string, string_array> fieldNames;
    fieldNames["node"].push_back( m_fieldName );

    CommunicationTools::
      SynchronizeFields( fieldNames,
      domain->getMeshBody( 0 )->getMeshLevel( 0 ),
      domain->getReference<array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
  }

  void LaplaceVEM::ApplyBoundaryConditions( real64 const time_n,
  real64 const dt,
  DomainPartition * const domain,
  DofManager const & dofManager,
  ParallelMatrix & matrix,
  ParallelVector & rhs )
  {
    ApplyDirichletBC_implicit( time_n + dt, dofManager, *domain, m_matrix, m_rhs );

    // Debug for logLevel >= 2
    GEOSX_LOG_LEVEL_RANK_0( 2, "After LaplaceVEM::ApplyBoundaryConditions" );
    GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
    GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );

    if( getLogLevel() >= 3 )
    {
      integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

      string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
      matrix.write( filename_mat, true );

      string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
      rhs.write( filename_rhs, true );

      GEOSX_LOG_RANK_0( "After LaplaceVEM::ApplyBoundaryConditions" );
      GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
      GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
    }
  }

  void LaplaceVEM::SolveSystem( DofManager const & dofManager,
  ParallelMatrix & matrix,
  ParallelVector & rhs,
  ParallelVector & solution )
  {
    rhs.scale( -1.0 ); // TODO decide if we want this here
    solution.zero();

    SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

    // Debug for logLevel >= 2
    GEOSX_LOG_LEVEL_RANK_0( 2, "After LaplaceVEM::SolveSystem" );
    GEOSX_LOG_LEVEL_RANK_0( 2, "\nSolution:\n" << solution );
  }

  void LaplaceVEM::ApplyDirichletBC_implicit( real64 const time,
  DofManager const & dofManager,
  DomainPartition & domain,
  ParallelMatrix & matrix,
  ParallelVector & rhs )
  {
    FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

    fsManager.Apply( time,
    &domain,
    "nodeManager",
    m_fieldName,
    [&]( FieldSpecificationBase const * const bc,
    string const &,
    set<localIndex> const & targetSet,
    Group * const targetGroup,
    string const GEOSX_UNUSED_ARG( fieldName ) )->void
    {
      bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( targetSet,
      time,
      targetGroup,
      m_fieldName,
      dofManager.getKey( m_fieldName ),
      1,
      matrix,
      rhs );
    });
  }
  //START_SPHINX_INCLUDE_00
  REGISTER_CATALOG_ENTRY( SolverBase, LaplaceVEM, std::string const &, Group * const )
  //END_SPHINX_INCLUDE_00
} /* namespace ANST */
